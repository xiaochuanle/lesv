#include "sv_signature.h"
#include "find_one_sv_group.h"
#include "sv_signature_file_name.h"
#include "sv_read_group_file_name.h"

#include "../../corelib/raw_reads_reader.h"

KHASH_SET_INIT_INT(int_set);

static const char* g_wrk_dir = NULL;
static const char* g_svsig_dir = NULL;
static const char* g_db_dir = NULL;

static void
s_dump_usage(const char* pn)
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s wrk_dir svsig_dir db_dir\n", pn);
}

static void
s_dump_outlier_reads(const char* db_dir,
    const int sid, 
    SvSignature* svsig_array, 
    const int svsig_count, 
    khash_t(int_set)* added_id_set)
{
    RawReadReader* raw_reads = RawReadReaderNew(db_dir, INIT_QUERY_DB_TITLE, TRUE);
    for (int i = 0; i < svsig_count; ++i) {
        int qid = svsig_array[i].qid;
        khiter_t pos = kh_get(int_set, added_id_set, qid);
        if (pos != kh_end(added_id_set)) continue;
        RawReadReader_SetReadFlag(raw_reads, qid);
    }
    RawReadReader_LoadFlaggedReads(raw_reads);

    char path[HBN_MAX_PATH_LEN];
    make_subject_outlier_fasta_path(g_wrk_dir, sid, path);
    hbn_dfopen(out, path, "w");
    kv_dinit(vec_u8, seq_v);
    ks_dinit(header);
    kh_clear(int_set, added_id_set);
    for (int i = 0; i < svsig_count; ++i) {
        int qid = svsig_array[i].qid;
        if (!raw_reads->raw_read_offset_array[qid]) continue;
        khiter_t pos = kh_get(int_set, added_id_set, qid);
        if (pos != kh_end(added_id_set)) continue;
        const char* name = RawReadReader_ReadName(raw_reads, qid);
        make_sv_read_header(name, svsig_array[i].qdir, sid, -1, svsig_array[i].fsfrom, svsig_array[i].fsto, &header);
        RawReadReader_ExtractRead(raw_reads, qid, FWD, &seq_v);
        u8* seq = kv_data(seq_v);
        int seq_l = kv_size(seq_v);
        for (int p = 0; p < seq_l; ++p) {
            int c = DECODE_RESIDUE(seq[p]);
            seq[p] = c;
        }
        fprintf(out, ">");
        hbn_fwrite(ks_s(header), 1, ks_size(header), out);
        fprintf(out, "\n");
        hbn_fwrite(seq, 1, seq_l, out);
        fprintf(out, "\n");
        int r = 0;
        kh_put(int_set, added_id_set, qid, &r);
        hbn_assert(r == 1);
    }

    ks_destroy(header);
    kv_destroy(seq_v);
    hbn_fclose(out);
    RawReadReaderFree(raw_reads);
}

static void
s_dump_ungrouped_svsig(SvSignature* sga, int sgc, void* added_id_set_)
{
    khash_t(int_set)* added_id_set = (khash_t(int_set)*)(added_id_set_);
    for (int i = 0; i < sgc; ++i) {
        khiter_t pos = kh_get(int_set, added_id_set, sga[i].qid);
        if (pos == kh_end(added_id_set)) {
            fprintf(stderr, "%d\t", i);
            dump_svsig(fprintf, stderr, sga[i], NULL);
        }
    }
}

static void
s_make_sv_read_group_for_one_subject(const char* db_dir, const int sid)
{
    if (subject_sv_read_group_is_created(g_wrk_dir, sid)) {
        HBN_LOG("sv signatures for subject %d are created. skip it.", sid);
        return;
    }

    char path[HBN_MAX_PATH_LEN];
    make_subject_sv_signature_path(g_svsig_dir, sid, path);
    int svsig_count = 0;
    SvSignature* svsig_array = NULL;
    int next_svsig_i = 0;
    int group_id = 0;
    char read_group_path[HBN_MAX_PATH_LEN];
    make_subject_sv_read_group_path(g_wrk_dir, sid, read_group_path);
    hbn_dfopen(out, read_group_path, "w");
    khash_t(int_set)* added_id_set = kh_init(int_set);

    svsig_array = load_sv_signature_array(path, &svsig_count, eGapAlignIns);
    HBN_LOG("load %d ins signatures", svsig_count);
    next_svsig_i = 0;
    while (next_svsig_i < svsig_count) {
        find_next_ins_group(svsig_array, svsig_count, &next_svsig_i, &group_id, out, added_id_set);
    }
    //s_dump_ungrouped_svsig(svsig_array, svsig_count, added_id_set);
    {
        next_svsig_i = 0;
        for (int i = 0; i < svsig_count; ++i) {
            if (svsig_array[i].type != eGapAlignInvalid) {
                svsig_array[next_svsig_i++] = svsig_array[i];
            }
        }
        svsig_count = next_svsig_i;
    }
    next_svsig_i = 0;
    while (next_svsig_i < svsig_count) {
        find_next_ins_group_relax(svsig_array, svsig_count, &next_svsig_i, &group_id, out, added_id_set);
    }
    free(svsig_array);
    svsig_count = 0;

    svsig_array = load_sv_signature_array(path, &svsig_count, eGapAlignDel);
    HBN_LOG("load %d del signatures", svsig_count);
    next_svsig_i = 0;
    while (next_svsig_i < svsig_count) {
        find_next_del_group(svsig_array, svsig_count, &next_svsig_i, &group_id, out, added_id_set);
    }
    {
        next_svsig_i = 0;
        for (int i = 0; i < svsig_count; ++i) {
            if (svsig_array[i].type != eGapAlignInvalid) {
                svsig_array[next_svsig_i++] = svsig_array[i];
            }
        }
        svsig_count = next_svsig_i;
    }
    next_svsig_i = 0;
    while (next_svsig_i < svsig_count) {
        find_next_del_group_relax(svsig_array, svsig_count, &next_svsig_i, &group_id, out, added_id_set);
    }
    //s_dump_ungrouped_svsig(svsig_array, svsig_count, added_id_set);
    free(svsig_array);
    svsig_count = 0;    

#if 0
    svsig_array = load_sv_signature_array(path, &svsig_count, eGapAlignInvalid);
    HBN_LOG("load %d total signatures", svsig_count);
    s_dump_outlier_reads(db_dir, sid, svsig_array, svsig_count, added_id_set);
    free(svsig_array);
    svsig_count = 0;
#endif

    kh_destroy(int_set, added_id_set);
    hbn_fclose(out);
    subject_sv_read_group_make_created(g_wrk_dir, sid);
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        s_dump_usage(argv[0]);
        return 1;
    }
    g_wrk_dir = argv[1];
    g_svsig_dir = argv[2];
    g_db_dir = argv[3];
    create_sv_read_group_dir(g_wrk_dir);
    const int num_subjects = seqdb_load_num_reads(g_db_dir, INIT_SUBJECT_DB_TITLE);
    for (int i = 0; i < num_subjects; ++i) {
        s_make_sv_read_group_for_one_subject(g_db_dir, i);
    }
}   
