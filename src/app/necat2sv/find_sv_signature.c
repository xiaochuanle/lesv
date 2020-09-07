#include "sv_reads.h"

#include "../../algo/hbn_traceback.h"
#include "../../algo/hbn_traceback_aux.h"
#include "../../algo/refine_align.h"
#include "../../corelib/line_reader.h"
#include "../../corelib/khash.h"
#include "../../corelib/m4_record.h"
#include "../../corelib/fasta.h"
#include "../../corelib/seq_tag_report.h"
#include "../../corelib/string2hsp.h"
#include "../../corelib/ksort.h"
#include "../../corelib/cstr_util.h"
#include "../../corelib/seqdb.h"
#include "../../corelib/raw_reads_reader.h"
#include "../../ncbi_blast/setup/hsp2string.h"
#include "../../corelib/cstr_util.h"
#include "../../ncbi_blast/setup/gapinfo.h"

#include "sv_signature.h"
#include "trf_array.h"
#include "sv_read_file_name.h"
#include "sv_signature_file_name.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <ctype.h>
#include <pthread.h>

static const char* g_wrk_dir = NULL;
static const char* g_svr_dir = NULL;
static const char* g_db_dir = NULL;
static int num_threads = -1;
static vec_sv_sig* ins_svsig_list = NULL;
static vec_sv_sig* del_svsig_list = NULL;
static pthread_mutex_t svsig_lock;
static SvRead* svr_array = NULL;
static int svr_count = 0;
static int svr_index = -1;
static pthread_mutex_t svr_index_lock;
static RawReadReader* query_reader = NULL;
static RawReadReader* subject_reader = NULL;
static int subject_id = -1;
static TrfArray* trf_array = NULL;
static int min_indel_size = 40;

static void
init_global_data(int argc, char* argv[])
{
    hbn_assert(argc == 6 || argc == 7);
    g_wrk_dir = argv[1];
    g_svr_dir = argv[2];
    g_db_dir = argv[3];
    min_indel_size = atoi(argv[4]);
    num_threads = atoi(argv[5]);
    const char* trf_array_path = (argc > 6) ? argv[6] : NULL;
    create_sv_signature_dir(g_wrk_dir);
    pthread_mutex_init(&svsig_lock, NULL);
    if (!ins_svsig_list) ins_svsig_list = (vec_sv_sig*)calloc(1, sizeof(vec_sv_sig));
    kv_init(*ins_svsig_list);
    if (!del_svsig_list) del_svsig_list = (vec_sv_sig*)calloc(1, sizeof(vec_sv_sig));
    kv_init(*del_svsig_list);
    pthread_mutex_init(&svr_index_lock, NULL);

    subject_reader = RawReadReaderNew(g_db_dir, INIT_SUBJECT_DB_TITLE, FALSE);
    trf_array = trf_array_path ? TrfArrayBuild(trf_array_path, g_db_dir, subject_reader->raw_read_info_array) : NULL;
}

static void
destroy_global_data()
{
    if (ins_svsig_list) {
        kv_destroy(*ins_svsig_list);
        free(ins_svsig_list);
        ins_svsig_list = NULL;
    }
    if (del_svsig_list) {
        kv_destroy(*del_svsig_list);
        free(del_svsig_list);
        del_svsig_list = NULL;
    }
    if (svr_array) {
        free(svr_array);
        svr_array = NULL;
    }
    if (query_reader) query_reader = RawReadReaderFree(query_reader);
    if (subject_reader) subject_reader = RawReadReaderFree(subject_reader);
    if (trf_array) trf_array = TrfArrayFree(trf_array);
}

static void
reset_global_data(const int sid)
{
    subject_id = sid;
    kv_clear(*ins_svsig_list);
    kv_clear(*del_svsig_list);
    if (svr_array) sfree(svr_array);
    char path[HBN_MAX_PATH_LEN];
    make_subject_sv_read_path(g_svr_dir, sid, path);
    svr_array = load_sv_read_array(path, &svr_count);
    svr_index = 0;
}

static int
get_next_svr_index()
{
    int i = -1;
    pthread_mutex_lock(&svr_index_lock);
    i = svr_index;
    ++svr_index;
    pthread_mutex_unlock(&svr_index_lock);
    if (i >= svr_count) i = -1;
    return i;
}

static void
add_svsig(vec_sv_sig* svsig_list, SvSignature* svsig)
{
    pthread_mutex_lock(&svsig_lock);
    kv_push(SvSignature, *svsig_list, *svsig);
    pthread_mutex_unlock(&svsig_lock);
}

void
find_sv_signature(
    const int query_id,
    const int qdir,
    const char* query_name,
    const char* qas,
    const char* sas,
    const int aln_size,
    const int subject_offset,
    int qb,
    int qe,
    int fqb,
    int fqe,
    int ql,
    int sb,
    int se,
    int fsb,
    int fse,
    int sl)
{
    const int E = min_indel_size;
    int i = 0;
    int qi = qb;
    int si = subject_offset + sb;
    SvSignature svpos;
    svpos.qid = query_id;
    svpos.qdir = qdir;
	static int ins_cnt = 0;
    while (i < aln_size) {
        if (qas[i] != GAP_CHAR && sas[i] != GAP_CHAR) {
            ++i;
            ++qi;
            ++si;
            continue;
        }

        if (qas[i] == GAP_CHAR) {
            hbn_assert(sas[i] != GAP_CHAR);
            int j = i + 1;
            while (j < aln_size && qas[j] == GAP_CHAR) {
                hbn_assert(sas[j] != GAP_CHAR);
                ++j;
            }
            int n = j - i;
            if (n >= E) {
                //HBN_LOG("find DEL pos at [%d, %d), len = %d", si, si + n, n);
                svpos.qfrom = qi;
                svpos.qto = qi + 1;
                svpos.fqfrom = fqb;
                svpos.fqto = fqe;
                svpos.qsize = ql;
                svpos.sfrom = si;
                svpos.sto = si + n;
                svpos.fsfrom = fsb;
                svpos.fsto = fse;
                svpos.ssize = sl;
                svpos.type = eGapAlignDel;
                add_svsig(del_svsig_list, &svpos);
            }
            si += n;
            i = j;
            continue;
        }

        if (sas[i] == GAP_CHAR) {
            hbn_assert(qas[i] != GAP_CHAR);
            int j = i + 1;
            while (j < aln_size && sas[j] == GAP_CHAR) {
                hbn_assert(qas[j] != GAP_CHAR);
                ++j;
            }
            int n = j - i;
            if (n >= E) {
                //HBN_LOG("find INS pos at %d, length = %d", si, n);
                svpos.qfrom = qi;
                svpos.qto = qi + n;
                svpos.fqfrom = fqb;
                svpos.fqto = fqe;
                svpos.qsize = ql;
                svpos.sfrom = si;
                svpos.sto = si + 1;
                svpos.fsfrom = fsb;
                svpos.fsto = fse;
                svpos.ssize = sl;
                svpos.type = eGapAlignIns;
                add_svsig(ins_svsig_list, &svpos);
            }
            qi += n;
            i = j;
            continue;
        }
    }
    hbn_assert(qi == qe);
    hbn_assert(si == subject_offset + se);
}

static void
s_dump_svsig_list(vec_sv_sig* svsig_list, FILE* out)
{
    SvSignature* svsiga = kv_data(*svsig_list);
    int svsigc = kv_size(*svsig_list);
    sort_svsig_sfrom_lt(svsigc, svsiga);

    for (int i = 0; i < svsigc; ++i) {
        const char* query_name = RawReadReader_ReadName(query_reader, svsiga[i].qid);
        dump_svsig(fprintf, out, svsiga[i], query_name);
    } 
}

void* discover_func(void* params)
{
    HbnTracebackData* tbck_data = HbnTracebackDataNew();
    kv_dinit(vec_u8, query_v);
    kv_dinit(vec_u8, subject_v);
    ks_dinit(qaln);
    ks_dinit(saln);
    while (1) {
        int i = get_next_svr_index();
        if (i == -1) break;
        //HBN_LOG("process query %d / %d", i, svr_count);
        SvRead svr = svr_array[i];
        int qid = svr.query_id;
        const char* query_name = RawReadReader_ReadName(query_reader, qid);
        int qdir = svr.qdir;
        int qoff = svr.qoff;
        int qend = svr.qend;
        int qsize = RawReadReader_ReadSize(query_reader, qid);
        if (svr.qdir == REV) {
            int x = qsize - qend;
            int y = qsize - qoff;
            qoff = x;
            qend = y;
        }
        int soff = svr.soff;
        int send = svr.send;
        //if (!(soff <= 15811790 && send >= 15811810)) continue;
        RawReadReader_ExtractSubRead(query_reader, qid, qoff, qend, qdir, &query_v);
        const u8* query = kv_data(query_v);
        int query_length = kv_size(query_v);
        if (send <= soff) dump_sv_read(fprintf, stderr, svr, NULL);
        RawReadReader_ExtractSubRead(subject_reader, subject_id, soff, send, FWD, &subject_v);
        const u8* subject = kv_data(subject_v);
        int subject_length = send - soff;
        const char* subject_name = RawReadReader_ReadName(subject_reader, subject_id);

        int qb, qe, sb, se;
        double perc_identity;
        int r = nw_ksw2_extd2(tbck_data->ksw,
                    qid,
                    query,
                    0,
                    query_length,
                    query_length,
                    subject_id,
                    subject,
                    0,
                    subject_length,
                    subject_length,
                    0,
                    0,
                    svr.dist * 1.2,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity,
                    &qaln,
                    &saln);
        if (!r) {
            HBN_LOG("normal ksw fail");
            r = nw_ksw2_extd2(tbck_data->ksw,
                    qid,
                    query,
                    0,
                    query_length,
                    query_length,
                    subject_id,
                    subject,
                    0,
                    subject_length,
                    subject_length,
                    0,
                    0,
                    hbn_max(query_length, subject_length),
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &perc_identity,
                    &qaln,
                    &saln);
        }
        if (!r) {
            dump_sv_read(fprintf, stderr, svr, query_name);
            HBN_LOG("ksw fail");
            continue;
        }

        int aln_size = ks_size(qaln);
        const char* qas = ks_s(qaln);
        const char* qae = qas + aln_size;
        const char* sas = ks_s(saln);
        const char* sae = sas + aln_size;
        r = truncate_align_bad_ends(qas, sas, aln_size, &qb, &qe, &sb, &se, &qas, &qae, &sas, &sae);
        if (!r) {
            dump_sv_read(fprintf, stderr, svr, query_name);
            HBN_LOG("truncate fail");
            continue;
        }
        aln_size = qae - qas;

        int fqb = qb + svr.qoff;
        int fqe = qe + svr.qoff;
        int fsb = sb + svr.soff;
        int fse = se + svr.soff;
        if (trf_array && fall_in_trf_subseq(trf_array, svr.subject_id, fsb, fse)) {
            HBN_LOG("****************[%d:%s, %d, %d] x [%d:%s, %d, %d] fall in trf subseq",
                svr.query_id, query_name, fqb, fqe, svr.subject_id, subject_name, fsb, fse);
            continue;
        }
        perc_identity = calc_ident_perc(qas, sas, aln_size, NULL, NULL);
        double eff_perc_identity = calc_effective_ident_perc(qas, sas, aln_size);
	    if (eff_perc_identity < 70.0) {
            HBN_LOG("xxxxxxxxxxxxxxxx[%d:%s, %d, %d] x [%d:%s, %d, %d] eff_perc_identity too low: %g, %g",
                svr.query_id, query_name, fqb, fqe, svr.subject_id, subject_name, fsb, fse, perc_identity, eff_perc_identity);
            continue;
        }
        query_length = RawReadReader_ReadSize(query_reader, svr.query_id);
        subject_length = RawReadReader_ReadSize(subject_reader, svr.subject_id);
        find_sv_signature(svr.query_id,
            svr.qdir,
            query_name,
            qas,
            sas,
            aln_size,
            soff,
            qb,
            qe,
            fqb,
            fqe,
            query_length,
            sb,
            se,
            fsb,
            fse,
            subject_length);        
    }        

    tbck_data = HbnTracebackDataFree(tbck_data);
    kv_destroy(query_v);
    kv_destroy(subject_v);
    ks_destroy(qaln);
    ks_destroy(saln);

    return NULL;
}

void
discover_one_subject(const int sid)
{
    if (subject_sv_signature_is_created(g_wrk_dir, sid)) {
        HBN_LOG("sv signature for subject %d is created, skip it.", sid);
        return;
    }
    reset_global_data(sid);

    query_reader = RawReadReaderNew(g_db_dir, INIT_QUERY_DB_TITLE, TRUE);
    for (int i = 0; i < svr_count; ++i) RawReadReader_SetReadFlag(query_reader, svr_array[i].query_id);
    RawReadReader_LoadFlaggedReads(query_reader);

    pthread_t job_ids[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(job_ids + i, NULL, discover_func, NULL);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(job_ids[i], NULL);
    }

    char path[HBN_MAX_PATH_LEN];
    make_subject_sv_signature_path(g_wrk_dir, sid, path);
    hbn_dfopen(out, path, "w");
    s_dump_svsig_list(ins_svsig_list, out);
    s_dump_svsig_list(del_svsig_list, out);
    hbn_fclose(out);
    query_reader = RawReadReaderFree(query_reader);
    subject_sv_signature_make_created(g_wrk_dir, sid);
}

static void
print_usage(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s wrk_dir svr_dir db_dir min_indel_size num_threads [trf_array_path]\n", pn);
}

int main(int argc, char* argv[])
{
    if (argc != 6 && argc != 7) {
        print_usage(argv[0]);
        return 1;
    }
    init_global_data(argc, argv);

    const int num_subjects = seqdb_load_num_reads(g_db_dir, INIT_SUBJECT_DB_TITLE);
    char job_name[256];
    char sid_str[64];
    for (int i = 0; i < num_subjects; ++i) {
        //if (i != 20) continue;
        u64_to_fixed_width_string_r(i, sid_str, HBN_DIGIT_WIDTH);
        sprintf(job_name, "detect signatures for subject %s", sid_str);
        hbn_timing_begin(job_name);
        discover_one_subject(i);
        hbn_timing_end(job_name);
    }
    destroy_global_data();

    return 0;
}
