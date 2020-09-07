#include "../../corelib/raw_reads_reader.h"
#include "../../algo/hbn_traceback.h"
#include "../../algo/hbn_traceback_aux.h"
#include "../necat2sv/sv_read_group_file_name.h"
#include "../../corelib/build_db.h"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/cns_read_header.h"
#include "map_results.h"
#include "sv_read_group_file_name.h"

#include <pthread.h>

typedef struct {
    RawReadReader* raw_reads;
    RawReadReader* db;
    int read_id_from;
    int read_id_to;
    int* read_id;
    pthread_mutex_t* read_id_lock;
    FILE* sam_out;
    pthread_mutex_t* out_lock;
    HbnTracebackData* tbck_data;
    int subject_id;
} ThreadWorkData;

ThreadWorkData*
ThreadWorkDataNew(RawReadReader* raw_reads,
    RawReadReader* db,
    int* read_id,
    pthread_mutex_t* read_id_lock,
    FILE* sam_out,
    pthread_mutex_t* out_lock,
    int subject_id)
{
    ThreadWorkData* data = (ThreadWorkData*)calloc(1, sizeof(ThreadWorkData));
    data->raw_reads = raw_reads;
    data->db = db;
    data->read_id_from = 1;
    data->read_id_to = 0;
    data->read_id = read_id;
    data->read_id_lock = read_id_lock;
    data->sam_out = sam_out;
    data->out_lock = out_lock;
    data->tbck_data = HbnTracebackDataNew();
    data->subject_id = subject_id;
    return data;
}

ThreadWorkData*
ThreadWorkDataFree(ThreadWorkData* data)
{
    HbnTracebackDataFree(data->tbck_data);
    free(data);
    return NULL;
}

static void
s_dump_sv_read_info(const int read_id,
    const int read_dir,
    const int read_cns_subseq_from,
    const int read_cns_subseq_to,
    const int read_length,
    const char* read_name,
    const u8* read,
    const char* subject_name,
    const int subject_id,
    const u8* subject,
    const int subject_subseq_offset,
    const int subject_total_ssize,
    const int qb,
    const int qe,
    const int sb,
    const int se,
    const char* qaln,
    const char* saln,
    const int aln_size,
    FILE* sam_out,
    pthread_mutex_t* out_lock)
{
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        read_id, read, qb, qe, qaln,
        subject_id, subject, sb, se, saln, aln_size, TRUE);
    //HBN_LOG("[%d, %d, %d, %d, %d] x [%d, %d]", read_lid, svri->fsqdir, qb, qe, svri->size, sb, se);
    int cns_qb = qb;
    int cns_qe = qe;
    if (read_dir == FWD) {
        cns_qb = read_cns_subseq_from;
        cns_qe = read_cns_subseq_to;
    } else {
        cns_qb = read_length - read_cns_subseq_to;
        cns_qe = read_length - read_cns_subseq_from;
    }

    int aln_from = 0;
    int qi = qb;
    if (cns_qb > qi) {
        while (qi < cns_qb) {
            if (qaln[aln_from] != GAP_CHAR) ++qi;
            ++aln_from;
        }
    }
    int aln_to = aln_size;
    qi = qe;
    if (qi > cns_qe) {
        while (qi > cns_qe) {
            if (qaln[aln_to - 1] != GAP_CHAR) --qi;
            --aln_to;
        }
    }
    if (aln_from >= aln_to) return;

    int qif = qb, qie = 0;
    int sif = sb, sie = 0;
    int ai = 0;
    while (ai < aln_from) {
        if (qaln[ai] != GAP_CHAR) ++qif;
        if (saln[ai] != GAP_CHAR) ++sif;
        ++ai;
    }
    hbn_assert(ai == aln_from);

    qie = qif;
    sie = sif;
    while (ai < aln_to) {
        if (qaln[ai] != GAP_CHAR) ++qie;
        if (saln[ai] != GAP_CHAR) ++sie;
        ++ai;
    }
    hbn_assert(ai == aln_to);

    //HBN_LOG("qif = %d, qie = %d, sif = %d, sie = %d, aln_from = %d, aln_to = %d",
    //    qif, qie, sif, sie, aln_from, aln_to);

    const char* q = qaln + aln_from;
    const char* s = saln + aln_from;
    const int as = aln_to - aln_from;
    //HBN_LOG("aln_from = %d, aln_to = %d, as = %d", aln_from, aln_to, as);
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        read_id, read, qif, qie, q,
        subject_id, subject, sif, sie, s, as, TRUE);

    int score = 0;
    double perc_identity = calc_ident_perc(q, s, as, NULL, &score);
    double eff_perc_identity = calc_effective_ident_perc(q, s, as);
    if (eff_perc_identity < 85.0) return;

#if 0
    HBN_LOG("[%d:%s, %d, %d, %d] x [%d, %d], perc_identity = %g, eff_perc_identity = %g", 
        read_id, read_name, qif, qie, read_length,
        sif, sie, 
        perc_identity, eff_perc_identity);
#endif

    ks_dinit(sam);
    ks_dinit(rg_info);
    make_sam_rg_info(subject_id, &rg_info);
    print_one_sam_result(read_name, read_dir, 0, qie - qif, qie - qif,
        subject_name, subject_subseq_offset + sif, subject_subseq_offset + sie, subject_total_ssize,
        q, s, as, score, perc_identity, eff_perc_identity, TRUE, ks_s(rg_info), &sam);

    pthread_mutex_lock(out_lock);
    if (sam_out) hbn_fwrite(ks_s(sam), 1, ks_size(sam), sam_out);
    pthread_mutex_unlock(out_lock);

    ks_destroy(rg_info);
    ks_destroy(sam);
}

static void*
map_sv_read_func(void* params)
{
    ThreadWorkData* data = (ThreadWorkData*)(params);
    const int subject_id = data->subject_id;
    const char* subject_name = RawReadReader_ReadName(data->db, subject_id);
    kv_dinit(vec_u8, read_v);
    kv_dinit(vec_u8, subject_v);

    while (1) {
        int read_id = -1;
        pthread_mutex_lock(data->read_id_lock);
        read_id = *data->read_id;
        ++(*data->read_id);
        pthread_mutex_unlock(data->read_id_lock);
        if (read_id >= data->read_id_to) break;

        const char* read_name = RawReadReader_ReadName(data->raw_reads, read_id);
        int read_dir, sid, gid, sfrom, sto;
        extract_info_from_sv_read_header(read_name, &read_dir, &sid, &gid, &sfrom, &sto);
        int read_cns_subseq_from, read_cns_subseq_to;
        extract_info_from_cns_read_header(read_name,
            NULL, NULL, &read_cns_subseq_from, &read_cns_subseq_to, NULL, NULL);
        //if (strcmp(read_name, 
        //    "02a176c7-3435-4f5d-b05f-70c52598b351_svr:1:20:935:10751024:10773136_cns:0:0:9282:11936:21505:21482")) {
        //    continue;
        //}
        //HBN_LOG("%s", read_name);
        //HBN_LOG("read_dir = %d, sid = %d, sfrom = %d, sto = %d, cns_from = %d, cns_to = %d",
        //    read_dir, sid, sfrom, sto, read_cns_subseq_from, read_cns_subseq_to);
        hbn_assert(sid == subject_id);
        const int read_length = RawReadReader_ReadSize(data->raw_reads, read_id);
        RawReadReader_ExtractRead(data->raw_reads, read_id, read_dir, &read_v);
        const u8* read = kv_data(read_v);
        RawReadReader_ExtractSubRead(data->db, subject_id, sfrom, sto, FWD, &subject_v);
        const u8* subject = kv_data(subject_v);
        const int subject_length = kv_size(subject_v);

        int qb, qe, sb, se;
        double perc_identity;
        int distance = hbn_max(read_length, subject_length);
        distance = distance * 0.2;
        kstring_t* qaln = &data->tbck_data->qabuf;
        kstring_t* saln = &data->tbck_data->sabuf;
        int r = nw_ksw2_extd2(data->tbck_data->ksw,
                    read_id,
                    read,
                    0,
                    read_length,
                    read_length,
                    subject_id,
                    subject,
                    0,
                    subject_length,
                    subject_length,
                    0,
                    0,
                    distance,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity,
                    qaln,
                    saln);
        if (!r) {
            distance = hbn_max(read_length, subject_length);
            r = nw_ksw2_extd2(data->tbck_data->ksw,
                    read_id,
                    read,
                    0,
                    read_length,
                    read_length,
                    subject_id,
                    subject,
                    0,
                    subject_length,
                    subject_length,
                    0,
                    0,
                    distance,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity,
                    qaln,
                    saln);            
        }
        if (!r) {
            HBN_LOG("ksw fail for %d", read_id);
            continue;
        }
        //dump_align_string(ks_s(*qaln), ks_s(*saln), ks_size(*qaln), stderr);
        ////HBN_LOG("[%d, %d, %d] x [%d, %d, %d], %g", qb, qe, read_length,
        //    sb, se, subject_length, perc_identity);

        s_dump_sv_read_info(read_id,
            read_dir,
            read_cns_subseq_from,
            read_cns_subseq_to,
            read_length,
            read_name,
            read,
            subject_name,
            subject_id,
            subject,
            sfrom,
            RawReadReader_ReadSize(data->db, subject_id),
            qb,
            qe,
            sb,
            se,
            ks_s(*qaln),
            ks_s(*saln),
            ks_size(*qaln),
            data->sam_out,
            data->out_lock);
        //exit(0);
    }
    kv_destroy(subject_v);
    kv_destroy(read_v);
    return NULL;
}

static void
s_make_subject_svr_path(const char* svr_raw_reads_dir, const int sid, char path[])
{
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    sprintf(path, "%s/s%s.svr.fasta", svr_raw_reads_dir, sid_str);
}

static void
s_make_subject_svr_pdb_path(const char* svr_raw_reads_dir, const int sid, char path[])
{
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    sprintf(path, "%s/s%s_svr_pdb", svr_raw_reads_dir, sid_str);
}

static RawReadReader*
s_make_subject_raw_reads(const char* svr_raw_reads_dir, const int subject_id)
{
    char svr_path[HBN_MAX_PATH_LEN];
    make_subject_corrected_normal_sv_read_path(svr_raw_reads_dir, subject_id, svr_path);
    char svr_pdb_path[HBN_MAX_PATH_LEN];
    s_make_subject_svr_pdb_path(svr_raw_reads_dir, subject_id, svr_pdb_path);
    create_directory(svr_pdb_path);
    build_db(svr_path, svr_pdb_path, INIT_QUERY_DB_TITLE, 0, INT32_MAX, 4000000000, 0);
    RawReadReader* raw_reads = RawReadReaderNew(svr_pdb_path, INIT_QUERY_DB_TITLE, FALSE);
    return raw_reads;
}

static void
s_map_sv_read_for_one_subject(RawReadReader* db, 
    const char* svr_raw_reads_dir,
    const int subject_id, 
    int num_threads, 
    FILE* sam_out)
{
    RawReadReader* raw_reads = s_make_subject_raw_reads(svr_raw_reads_dir, subject_id);
    ThreadWorkData* data[num_threads];
    int read_id = 0;
    pthread_mutex_t read_id_lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;
    for (int i = 0; i < num_threads; ++i) {
        data[i] = ThreadWorkDataNew(raw_reads,
                    db,
                    &read_id,
                    &read_id_lock,
                    sam_out,
                    &out_lock,
                    subject_id);
    }
    int next_batch_id = 0;
    const int kBatchSize = 5000;
    pthread_t jobs[num_threads];
    while (1) {
        int from = next_batch_id;
        if (from >= raw_reads->dbinfo.num_seqs) break;
        int to = hbn_min(from + kBatchSize, raw_reads->dbinfo.num_seqs);
        read_id = from;
        for (int i = 0; i < num_threads; ++i) {
            data[i]->read_id_from = from;
            data[i]->read_id_to = to;
            pthread_create(jobs + i, NULL, map_sv_read_func, data[i]);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(jobs[i], NULL);
        }
        next_batch_id = to;
    }
    for (int i = 0; i < num_threads; ++i) {
        data[i] = ThreadWorkDataFree(data[i]);
    }
    raw_reads = RawReadReaderFree(raw_reads);
}

static void
print_usage(const char* pn)
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s subject_pdb_dir svr_group_dir num_threads\n", pn);
}

static FILE*
s_open_sam_file_for_mapping_subject_sv_reads(const char* svr_raw_reads_dir, 
    int subject_id,
    RawReadReader* db,
    int argc,
    char* argv[])
{
    char path[HBN_MAX_PATH_LEN];
    make_subject_corrected_normal_sv_read_path(svr_raw_reads_dir, subject_id, path);
    strcat(path, ".sam");
    hbn_dfopen(out, path, "w");

    ks_dinit(rg_info);
    make_sam_rg_info(subject_id, &rg_info);
    print_sam_prolog(out, kSamVersion, HBN_PACKAGE_VERSION, ks_s(rg_info), 
        DEFAULT_SAMPLE, NULL, db, subject_id, argc, argv);
    ks_destroy(rg_info);

    return out;
}

static BOOL
s_subject_sv_reads_are_mapped(const char* svr_raw_reads_dir, int subject_id)
{
    char path[HBN_MAX_PATH_LEN];
    make_subject_corrected_normal_sv_read_path(svr_raw_reads_dir, subject_id, path);
    strcat(path, ".mapped");
    return access(path, F_OK) == 0; 
}

static void
s_subject_sv_reads_make_mapped(const char* svr_raw_reads_dir, int subject_id)
{
    char path[HBN_MAX_PATH_LEN];
    make_subject_corrected_normal_sv_read_path(svr_raw_reads_dir, subject_id, path);
    strcat(path, ".mapped");
    hbn_dfopen(out, path, "w");
    hbn_fclose(out);
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }
    const char* subject_pdb_dir = argv[1];
    const char* svr_raw_reads_dir = argv[2];
    const int num_threads = atoi(argv[3]);
    RawReadReader* db = RawReadReaderNew(subject_pdb_dir, INIT_SUBJECT_DB_TITLE, FALSE);
    char job_name[256];
    for (int i = 0; i < db->dbinfo.num_seqs; ++i) {
        //if (i != 20) continue;
        char sid_str[64];
        u64_to_fixed_width_string_r(i, sid_str, HBN_DIGIT_WIDTH);
        if (s_subject_sv_reads_are_mapped(svr_raw_reads_dir, i)) {
            HBN_LOG("sv reads are subject %s have been mapped.", sid_str);
            continue;
        }
        sprintf(job_name, "mapping sv reads for subject %s", sid_str);
        hbn_timing_begin(job_name);
        FILE* sam_out = s_open_sam_file_for_mapping_subject_sv_reads(svr_raw_reads_dir, 
                            i,
                            db,
                            argc,
                            argv);
        s_map_sv_read_for_one_subject(db, svr_raw_reads_dir, i, num_threads, sam_out);
        hbn_fclose(sam_out);
        hbn_timing_end(job_name);
        s_subject_sv_reads_make_mapped(svr_raw_reads_dir, i);
    }
    RawReadReaderFree(db);
    return 0;
}