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
#include "../../corelib/partition_mt.h"
#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/setup/hsp2string.h"
#include "../../corelib/cstr_util.h"
#include "../../corelib/raw_reads_reader.h"

#include "align_subseqs.h"
#include "trf_array.h"
#include "sv_read_file_name.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <ctype.h>
#include <pthread.h>

pthread_mutex_t gq_lock;

static int g_min_seq_size = 3000;
static double g_min_ident_perc = 70.0;
static int g_max_overhang = 300;
static const char* g_wrk_dir = NULL;
static const char* g_db_dir = NULL;

static void
set_sv_read_from_m4(const M4Record* m4, SvRead* rp)
{
    rp->query_id = m4->qid;
    rp->qdir = m4->qdir;
    rp->qsize = m4->qsize;
    rp->qoff = m4->qoff;
    rp->qend = m4->qend;
    rp->subject_id = m4->sid;
    rp->soff = m4->soff;
    rp->send = m4->send;
    int dist = (m4->qend - m4->qoff + m4->send - m4->soff) * (100.0 - m4->ident_perc) / 200.0;
    rp->dist = dist;
    rp->dual_qdir = -1;
}

static void
s_dump_one_item_info(const char* title, int seq, size_t res)
{
    char buf1[64], buf2[64];
    u64_to_string_comma(seq, buf1);
    u64_to_string_datasize(res, buf2);
    fprintf(stderr, "%s: %s (%s)\n", title, buf1, buf2);  
}

typedef struct {
    int query_vol_index;
    RawReadReader* query_reader;
    RawReadReader* subject_reader;
    vec_sv_read indel_sv_read_list;
    HbnTracebackData* tbck_data;
    M4Record* m4a;
    int m4c;
    int* m4_idx;
    pthread_mutex_t* m4_idx_lock;
    SubseqAlignResult* align_result;
    TrfArray* trf_array;
} WorkData;

WorkData*
WorkDataNew()
{
    WorkData* data = (WorkData*)calloc(1, sizeof(WorkData));
    data->query_vol_index = -1;
    data->query_reader = NULL;
    data->subject_reader = NULL;
    kv_init(data->indel_sv_read_list);
    data->tbck_data = HbnTracebackDataNew();
    data->m4a = NULL;
    data->m4c = 0;
    data->m4_idx = NULL;
    data->m4_idx_lock = NULL;
    data->align_result = SubseqAlignResultNew();
    data->trf_array = NULL;
    return data;
}

WorkData*
WorkDataFree(WorkData* data)
{
    if (!data) return NULL;
    kv_destroy(data->indel_sv_read_list);
    data->tbck_data = (data->tbck_data) ? HbnTracebackDataFree(data->tbck_data) : NULL;
    data->align_result = SubseqAlignResultFree(data->align_result);
    sfree(data);
    return NULL;
}

static void
WorkData_ExtractQuerySubseq(WorkData* data, int query_gid, int from, int to, int dir, vec_u8* query)
{
    hbn_assert(data->query_vol_index >= 0);
    RawReadReader_ExtractSubRead(data->query_reader, query_gid, from, to, dir, query);
}

const char*
WorkData_ExtractQueryName(WorkData* data, int query_gid)
{
    hbn_assert(data->query_vol_index >= 0);
    return RawReadReader_ReadName(data->query_reader, query_gid);
}

static BOOL
m4_is_complete(M4Record* m4, int max_overhang)
{
    int qoff = m4->qoff;
    int qend = m4->qend;
    int qsize = m4->qsize;
    return (qoff <= max_overhang) && (qsize - qend <= max_overhang);    
}

static BOOL
calc_m4_eff_perc_identity(WorkData* data, M4Record* m4, double* eff_perc_identity, int* dist, int* score)
{
    int qb = m4->qoff;
    int qe = m4->qend;
    if (m4->qdir == REV) {
        int x = m4->qsize - m4->qend;
        int y = m4->qsize - m4->qoff;
        qb = x;
        qe = y;
    }
    int sb = m4->soff;
    int se = m4->send;
    BOOL r = FALSE;
    int max_dist = hbn_max(m4->qend - m4->qoff, m4->send - m4->soff);
    max_dist = max_dist * (100.0 - m4->ident_perc) / 100.0 * 1.2;
    //r = align_and_refine_subseq_with_edlib(data->query_reader, data->subject_reader, data->can_finder,
    //        data->tbck_data, max_dist, m4->qid, m4->qdir, qb, qe, m4->sid, sb, se, data->align_result);
    if (!r) {
        r = align_and_refine_subseq_with_ksw(data->query_reader, data->subject_reader,
            data->tbck_data, max_dist, m4->qid, m4->qdir, qb, qe, m4->sid, sb, se, data->align_result);
        //if (r) HBN_LOG("rescue, [%d, %d, %d] x [%d, %d, %d], %g, %g", m4->qid, data->align_result->qb, data->align_result->qe,
        //            m4->sid, data->align_result->sb, data->align_result->se, data->align_result->perc_identity, data->align_result->eff_perc_identity);
        if (!r) {
            DUMP_M4_RECORD(fprintf, stderr, *m4);
            HBN_LOG("ksw fail, max_dist = %d", max_dist);
        }
    }
    if (r) {
        if (eff_perc_identity) *eff_perc_identity = data->align_result->eff_perc_identity;
        if (dist) *dist = data->align_result->dist;
        if (score) *score = data->align_result->score;
    }
    return r;
}

#define m4_qcov_gt(a, b) ((a).qend - (a).qoff > (b).qend - (b).qoff)
KSORT_INIT(m4_qcov_gt, M4Record, m4_qcov_gt);

static BOOL
find_complete_m4(WorkData* data, M4Record* m4a, int m4c)
{
    M4Record cmp_m4a[m4c];
    int cmp_cnt = 0;
    for (int i = 0; i < m4c; ++i) if (m4_is_complete(m4a + i, g_max_overhang)) cmp_m4a[cmp_cnt++] = m4a[i];
    if (!cmp_cnt) return FALSE;

    if (cmp_cnt == 1 && cmp_m4a[0].ident_perc >= g_min_ident_perc) {
        M4Record* m4 = cmp_m4a;
        if (data->trf_array 
            &&
            fall_in_trf_subseq(data->trf_array, m4->sid, m4->soff, m4->send)) {
            //DUMP_M4_RECORD(fprintf, stderr, *m4);
            //HBN_LOG("fall in trf subseq");
            return TRUE;
        }
        SvRead rp;
        set_sv_read_from_m4(m4, &rp);
        kv_push(SvRead, data->indel_sv_read_list, rp);
        return TRUE;
    }

    for (int i = 0; i < cmp_cnt; ++i) {
        double eff_perc_identity = 0.0;
        int r = calc_m4_eff_perc_identity(data, m4a + i, &eff_perc_identity, NULL, NULL);
        if (!r) cmp_m4a[i].qid = -1;
        if (eff_perc_identity < g_min_ident_perc) cmp_m4a[i].qid = -1;
        cmp_m4a[i].ident_perc = eff_perc_identity;
    }

    int m = 0;
    for (int i = 0; i < cmp_cnt; ++i) {
        if (cmp_m4a[i].qid != -1) {
            //HBN_LOG("high perc identity m4 %d", i);
            //DUMP_M4_RECORD(fprintf, stderr, cmp_m4a[i]);
            cmp_m4a[m++] = cmp_m4a[i];
        }
    }
    cmp_cnt = m;
    if (!cmp_cnt) return TRUE;

    if (cmp_cnt == 1) {
        M4Record* m4 = cmp_m4a;
        if (data->trf_array 
            &&
            fall_in_trf_subseq(data->trf_array, m4->sid, m4->soff, m4->send)) {
            //DUMP_M4_RECORD(fprintf, stderr, *m4);
            //HBN_LOG("fall in trf subseq");
            return TRUE;
        }
        SvRead rp;
        set_sv_read_from_m4(m4, &rp);
        kv_push(SvRead, data->indel_sv_read_list, rp);
        //HBN_LOG("find unique complete m4");
        //DUMP_M4_RECORD(fprintf, stderr, *m4);
        return TRUE;
    }

    ks_introsort_m4_ident_gt(cmp_cnt, cmp_m4a);
    if (cmp_m4a[0].ident_perc - cmp_m4a[1].ident_perc > 10.0) {
        M4Record* m4 = cmp_m4a;
        if (data->trf_array 
            &&
            fall_in_trf_subseq(data->trf_array, m4->sid, m4->soff, m4->send)) {
            //DUMP_M4_RECORD(fprintf, stderr, *m4);
            //HBN_LOG("fall in trf subseq");
            return TRUE;
        }
        SvRead rp;
        set_sv_read_from_m4(m4, &rp);
        kv_push(SvRead, data->indel_sv_read_list, rp);
        //HBN_LOG("find best complete m4");
        //DUMP_M4_RECORD(fprintf, stderr, *m4);
        return TRUE;        
    }
    return TRUE;
}

static BOOL 
two_m4s_are_dual(M4Record* lm, M4Record* rm)
{
    int lsb, lse, rsb, rse;
    if (lm->soff < rm->soff) {
        lsb = lm->soff;
        lse = lm->send;
        rsb = rm->soff;
        rse = rm->send;
    } else {
        lsb = rm->soff;
        lse = rm->send;
        rsb = lm->soff;
        rse = lm->send;
    }
    hbn_assert(lsb <= rsb);

    if (rsb <= lse) return TRUE; // overlap
    hbn_assert(rsb > lse);
    int dist_s = rsb - lse;
    int r = (dist_s <= 30000);

    return r;
}

static void
find_dual_m4s_on_one_subject(WorkData* data, M4Record* m4a, int m4c, vec_m4* dual_m4_list)
{
    M4Record lm[m4c];
    int lc = 0;
    M4Record rm[m4c];
    int rc = 0;
    for (int i = 0; i < m4c; ++i) {
        M4Record* m4 = m4a + i;
        if (m4_is_complete(m4, g_max_overhang)) continue;
        if (data->trf_array 
            &&
            fall_in_trf_subseq(data->trf_array, m4->sid, m4->soff, m4->send)) {
            //DUMP_M4_RECORD(fprintf, stderr, *m4);
            //HBN_LOG("fall in trf subseq");
            continue;
        }
        if (m4->qoff <= g_max_overhang) lm[lc++] = *m4;
        if (m4->qsize - m4->qend <= g_max_overhang) rm[rc++] = *m4;
    }
    if (lc == 0 || rc == 0) return;

    for (int i = 0; i < lc; ++i) {
        M4Record* m4 = lm + i;
        if (m4->ident_perc < g_min_ident_perc) {
            double eff_perc_identity;
            int r = calc_m4_eff_perc_identity(data, m4, &eff_perc_identity, NULL, NULL);
            if (r == 0 || eff_perc_identity < g_min_ident_perc) m4->qid = -1;
            m4->ident_perc = eff_perc_identity;
        }
    }
    int m = 0;
    for (int i = 0; i < lc; ++i) {
        if (lm[i].qid != -1) lm[m++] = lm[i];
    }
    lc = m;
    if (!lc) return;

    for (int i = 0; i < rc; ++i) {
        M4Record* m4 = rm + i;
        if (m4->ident_perc < g_min_ident_perc) {
            double eff_perc_identity;
            int r = calc_m4_eff_perc_identity(data, m4, &eff_perc_identity, NULL, NULL);
            if (r == 0 || eff_perc_identity < g_min_ident_perc) m4->qid = -1;
            m4->ident_perc = eff_perc_identity;
        }
    }
    m = 0;
    for (int i = 0; i < rc; ++i) {
        if (rm[i].qid != -1) rm[m++] = rm[i];
    }
    rc = m;
    if (!rc) return;

    int nd = 0;
    for (int i = 0; i < lc; ++i) {
        M4Record* left = lm + i;
        for (int j = 0; j < rc; ++j) {
            M4Record* right = rm + j;
            if (two_m4s_are_dual(left, right)) {
                kv_push(M4Record, *dual_m4_list, *left);
                kv_push(M4Record, *dual_m4_list, *right);
                ++nd;
            }
        }
    }
}

static BOOL
s_chain_dual_m4s(WorkData* data, M4Record* m1, M4Record* m2)
{
//if (m1->qid != 18565) return;
    //HBN_LOG("chaining two m4");
    //DUMP_M4_RECORD(fprintf, stderr, *m1);
    //DUMP_M4_RECORD(fprintf, stderr, *m2);
    if (m1->qdir != m2->qdir) {
        //HBN_LOG("different qdir");
        return FALSE;
    }
    int lqb, lqe, lsb, lse;
    int rqb, rqe, rsb, rse;
    if (m1->qdir == FWD) {
            lqb = m1->qoff;
            lqe = m1->qend;
            lsb = m1->soff;
            lse = m1->send;
            
            rqb = m2->qoff;
            rqe = m2->qend;
            rsb = m2->soff;
            rse = m2->send;
    } else {
            lqb = m2->qsize - m2->qend;
            lqe = m2->qsize - m2->qoff;
            lsb = m2->soff;
            lse = m2->send;
            
            rqb = m1->qsize - m1->qend;
            rqe = m1->qsize - m1->qoff;
            rsb = m1->soff;
            rse = m1->send;
    }
    if (lsb > rsb) {
        //HBN_LOG("wrong soff");
        return FALSE;
    }   

    int qoff = hbn_min(lqb, rqb);
    int qend = hbn_max(lqe, rqe);
    int soff = hbn_min(lsb, rsb);
    int send = hbn_max(lse, rse);
    const char* query_name = WorkData_ExtractQueryName(data, m1->qid);
    int qb, qe, sb, se;
    double perc_identity;
    double eff_perc_identity;

    if (!align_and_refine_subseq_with_ksw(data->query_reader, data->subject_reader, 
        data->tbck_data, -1, m1->qid, m1->qdir, qoff, qend, m1->sid, 
        soff, send, data->align_result)) {
	    return FALSE;
    }
    qb = data->align_result->qb;
    qe = data->align_result->qe;
    sb = data->align_result->sb;
    se = data->align_result->se;
    perc_identity = data->align_result->perc_identity;
    eff_perc_identity = data->align_result->eff_perc_identity;
    
    //HBN_LOG("qid = %d, sid = %d, perc_identity = %g, eff_perc_identity = %g",m1->qid, m1->sid, perc_identity, eff_perc_identity);
    if (eff_perc_identity > m1->ident_perc - 4.0 || eff_perc_identity > m2->ident_perc - 4.0) {
#if 0
        pthread_mutex_lock(&gq_lock);
        HBN_LOG("extending %d:%s", m1->qid, query_name);
        DUMP_M4_RECORD(fprintf, stderr, *m1);
        DUMP_M4_RECORD(fprintf, stderr, *m2);
        HBN_LOG("[%d, %d, %d] x [%d, %d, %d], perc_identity = %g, eff_perc_identity = %g", qb, qe, m1->qsize, sb, se, m1->ssize, perc_identity, eff_perc_identity);
        pthread_mutex_unlock(&gq_lock); 
#endif
        SvRead rp; memset(&rp, 0, sizeof(SvRead));
        init_sv_read(rp);
        rp.query_id = m1->qid;
        rp.qdir = m1->qdir;
        if (rp.qdir == FWD) {
            rp.qoff = qb;
            rp.qend = qe;
        } else {
            rp.qoff = m1->qsize - qe;
            rp.qend = m1->qsize - qb;
        }
        rp.qsize = m1->qsize;
        rp.subject_id = m1->sid;
        rp.soff = sb;
        rp.send = se;
        rp.dist = data->align_result->dist;
        kv_push(SvRead, data->indel_sv_read_list, rp);
    } 

    return TRUE;
}

static BOOL 
find_dual_m4s(WorkData* data, M4Record* m4a, int m4c)
{
    kv_dinit(vec_m4, dual_m4_list);
    int i = 0;
    while (i < m4c) {
        int sid = m4a[i].sid;
        int j = i + 1;
        while (j < m4c && m4a[j].sid == sid) ++j;
        find_dual_m4s_on_one_subject(data, m4a + i, j - i, &dual_m4_list);
        i = j;
    }
    int nd = kv_size(dual_m4_list);
    hbn_assert((nd%2) == 0);
    nd /= 2;
    if (nd == 1) {
        M4Record* m1 = &kv_A(dual_m4_list, 0);
        M4Record* m2 = &kv_A(dual_m4_list, 1);
        hbn_assert(m1->sid == m2->sid);
	    BOOL chain_dual_m4_success = FALSE;
        chain_dual_m4_success = s_chain_dual_m4s(data, m1, m2);
    }
    kv_destroy(dual_m4_list);
    return (nd > 0);
}

static BOOL
find_valid_incomplete_m4(WorkData* data, M4Record* m4a, int m4c)
{
    if (m4c == 0) return FALSE;
    ks_introsort_m4_qcov_gt(m4c, m4a);
    if (m4c > 1) {
        M4Record* m1 = m4a;
        M4Record* m2 = m4a + 1;
        int len1 = m1->qend - m1->qoff;
        int len2 = m2->qend - m2->qoff;
        if (len2 > len1 * 0.8) return FALSE;
    }
    M4Record* m4 = m4a;
    int len = m4->qend - m4->qoff;
    if (len < m4->qsize * 0.8) return FALSE;

    if (m4->ident_perc < g_min_ident_perc) {
        double eff_perc_identity;
        int r = calc_m4_eff_perc_identity(data, m4, &eff_perc_identity, NULL, NULL);
        if (!r) return FALSE;
        if (eff_perc_identity < g_min_ident_perc) return FALSE;
        m4->ident_perc = eff_perc_identity;
    }

    SvRead rp;
    set_sv_read_from_m4(m4, &rp);
    kv_push(SvRead, data->indel_sv_read_list, rp);
    //HBN_LOG("find incomplete hsp");
    //DUMP_M4_RECORD(fprintf, stderr, *m4);
    return TRUE;
}

static void
s_remove_contained_m4s_one_subject(M4Record* m4a, int m4c)
{
    const int E = 200;
    for (int i = 0; i < m4c; ++i) {
        M4Record* mi = m4a + i;
        if (mi->qid == -1) continue;
        int qoff = mi->qoff;
        int qend = mi->qend;
        int soff = mi->soff;
        int send = mi->send;
        for (int j = i + 1; j < m4c; ++j) {
            M4Record* mj = m4a + j;
            if (mj->qid == -1) continue;
            if (mi->qdir != mj->qdir) continue;
            int qb = mj->qoff;
            int qe = mj->qend;
            int sb = mj->soff;
            int se = mj->send;
            int r = (qb + E >= qoff && qe <= qend + E) && (sb + E >= soff && se <= send + E);
            if (r) {
                mj->qid = -1;
                continue;
            }
            r = (qoff + E >= qb && qend <= qe + E) && (soff + E >= sb && send <= se + E);
            if (r) {
                mi->qid = -1;
                continue;
            }
        }
    }
}

void
remove_contained_m4s(M4Record* m4a, int* m4c_)
{
    int cnt = *m4c_;
    int i = 0;
    while (i < cnt) {
        int sid = m4a[i].sid;
        int j = i + 1;
        while (j < cnt) {
            int sid1 = m4a[j].sid;
            if (sid1 != sid) break;
            ++j;
        }
        s_remove_contained_m4s_one_subject(m4a + i, j - i);
        i = j;
    }
    int m = 0;
    for (i = 0; i < cnt; ++i) if (m4a[i].qid != -1) m4a[m++] = m4a[i];
    *m4c_ = m;
}

void
remove_repeat_m4s(M4Record* m4a, int* m4c_)
{
    const int E = 300;
    int m4c = *m4c_;
    for (int i = 0; i < m4c; ++i) {
        M4Record* mi = m4a + i;
        if (mi->qid == -1) continue;
        int qoff = mi->qoff;
        int qend = mi->qend;
        for (int j = i + 1; j < m4c; ++j) {
            M4Record* mj = m4a + j;
            if (mj->qid == -1) continue;
            int qb = mj->qoff;
            int qe = mj->qend;

            int max_qoff = hbn_max(qoff, qb);
            int min_qend = hbn_min(qend, qe);
            if (max_qoff < min_qend) {
                int x = (max_qoff >= qoff) ? (max_qoff - qoff) : (qoff - max_qoff);
                int y = (min_qend >= qend) ? (min_qend - qend) : (qend - min_qend);
                int u = (max_qoff >= qb) ? (max_qoff - qb) : (qb - max_qoff);
                int v = (min_qend >= qe) ? (min_qend - qe) : (qe - min_qend);
                if (x + y <= E && u + v <= E) {
#if 0
                    HBN_LOG("find repeat hsps:");
                    dump_blasthsp(fprintf, stderr, *hsp);
                    dump_blasthsp(fprintf, stderr, *hsp1);
#endif
                    mi->qid = -1;
                    mj->qid = -1;
                }
                continue;
            }
        }        
    }

    int m = 0;
    for (int i = 0; i < m4c; ++i) if (m4a[i].qid != -1) m4a[m++] = m4a[i];
    *m4c_ = m;
}

static void
s_find_sv_signature(
    const char* query_name,
    const char* qas,
    const char* sas,
    const int aln_size,
    const int subject_offset,
    int qb,
    int qe,
    int sb,
    int se)
{
    const int E = 40;
    int i = 0;
    int qi = qb;
    int si = subject_offset + sb;
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
                HBN_LOG("**** find DEL pos at [%d, %d), len = %d", si, si + n, n);
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
                HBN_LOG("**** find INS pos at %d, length = %d", si, n);
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
process_specific_m4_list(WorkData* data, M4Record* m4a, int m4c)
{
    const char* target1 = "80f02c09-1ecd-4568-a256-93bffd841e80";
    const char* target2 = "20e4c353-6b85-4e8e-86e3-63984260199d";
    const char* name = RawReadReader_ReadName(data->query_reader, m4a[0].qid);
    int r = (strcmp(target1, name) == 0) || (strcmp(target2, name) == 0);
    if (!r) return;
    static int cnt = 0;
    for (int i = 0; i < m4c; ++i) {
        M4Record* m4 = m4a + i;
        int qid = m4->qid;
        int qdir = m4->qdir;
        int qoff = m4->qoff;
        int qend = m4->qend;
        int qsize = m4->qsize;
        if (m4->qdir == REV) {
            int x = qsize - qend;
            int y = qsize - qoff;
            qoff = x;
            qend = y;
        }
        int subject_id = m4->sid;
        int soff = m4->soff;
        int send = m4->send;
        //if (!(m4->ssize == 48129895 && soff <= 19326600 && send >= 19329110)) continue;
	    DUMP_M4_RECORD(fprintf, stderr, *m4);
        int r = align_and_refine_subseq_with_ksw(data->query_reader,
                    data->subject_reader,
                    data->tbck_data,
                    -1,
                    qid,
                    qdir,
                    qoff,
                    qend,
                    subject_id,
                    soff,
                    send,
                    data->align_result);
        if (!r) continue;
        dump_align_string(data->align_result->qas, data->align_result->sas, data->align_result->qae - data->align_result->qas, stderr);
        const char* qname = RawReadReader_ReadName(data->query_reader, qid);
        const char* sname = RawReadReader_ReadName(data->subject_reader, subject_id);
        HBN_LOG("%d, [%d:%s, %d, %d, %d] x [%d, %s, %d, %d, %d], %g, %g", cnt, qid, qname, data->align_result->qb, data->align_result->qe, m4->qsize,
            subject_id, sname, data->align_result->sb, data->align_result->se, m4->ssize, data->align_result->perc_identity, data->align_result->eff_perc_identity);
        ++cnt;
        s_find_sv_signature(qname, 
            data->align_result->qas,
            data->align_result->sas,
            data->align_result->qae - data->align_result->qas,
            0,
            data->align_result->qb,
            data->align_result->qe,
            data->align_result->sb,
            data->align_result->se);
    }     
}

static void
process_one_m4_list(WorkData* data, M4Record* m4a, int m4c)
{
    //process_specific_m4_list(data, m4a, m4c);
    //return;
    int query_length = m4a[0].qsize;
    if (query_length < g_min_seq_size) return;
    int qid = m4a[0].qid;
    //const char* query_name = WorkData_ExtractQueryName(data, m4a[0].qid);
    //if (strcmp(query_name, "db6c8deb-9a3a-4db2-850c-a33ce8754bc1")) return;
    //HBN_LOG("processing query %d:%s", qid, query_name);
    remove_contained_m4s(m4a, &m4c);
    hbn_assert(m4c);
    //HBN_LOG("find complete_hsps, %d", m4c);
    if (find_complete_m4(data, m4a, m4c)) return;
    remove_repeat_m4s(m4a, &m4c);
    if (!m4c) return;
    //HBN_LOG("find dual hsps, %d", m4c);
    if (find_dual_m4s(data, m4a, m4c)) return;
    //HBN_LOG("find valid incomplete hsps");
    //if (find_valid_incomplete_m4(data, m4a, m4c)) return;
}

static BOOL 
get_next_m4a(WorkData* data, M4Record** m4ap, int* m4cp)
{
    int m4i = -1;
    int m4c = 0;
    pthread_mutex_lock(data->m4_idx_lock);
    m4i = *data->m4_idx;
    if (m4i < data->m4c) {
        int qid = data->m4a[m4i].qid;
        *m4ap = data->m4a + m4i;
        while (m4i < data->m4c) {
            if (data->m4a[m4i].qid != qid) break;
            ++m4i;
            ++m4c;
        }
    }
    *data->m4_idx = m4i;
    pthread_mutex_unlock(data->m4_idx_lock);
    *m4cp = m4c;
    return (m4c > 0);
}

static void*
thread_func(void* params)
{
    WorkData* data = (WorkData*)(params);
    M4Record* m4a;
    int m4c;
    while (get_next_m4a(data, &m4a, &m4c)) {
        process_one_m4_list(data, m4a, m4c);
    }
    return NULL;
}

static void
load_volume_mapping_results(const char* db_dir, int qi, int sj, M4Record** m4ap, size_t* m4cp)
{
    char path[HBN_MAX_PATH_LEN];
    char qi_buf[64];
    char sj_buf[64];
    u64_to_fixed_width_string_r(qi, qi_buf, HBN_DIGIT_WIDTH);
    u64_to_fixed_width_string_r(sj, sj_buf, HBN_DIGIT_WIDTH);
    sprintf(path, "%s/backup_results/stagebackup_results_Q%s_D%s.m4x", db_dir, qi_buf, sj_buf);
    *m4ap = (M4Record*)load_part_records(path, sizeof(M4Record), m4cp);
}

static void
dump_query_vol_sv_read_results(WorkData** wd_array, int num_threads, int qi)
{
    char path[HBN_MAX_PATH_LEN];
    make_query_vol_sv_read_path(g_wrk_dir, qi, path);
    hbn_dfopen(out, path, "w");
    for (int i = 0; i < num_threads; ++i) {
        WorkData* data = wd_array[i];
        for (size_t p = 0; p < kv_size(data->indel_sv_read_list); ++p) {
            SvRead rp = kv_A(data->indel_sv_read_list, p);
            dump_sv_read(fprintf, out, rp, NULL);
        }
    }
    hbn_fclose(out);
}

static RawReadReader*
s_load_query_reader(const char* db_dir, int vid)
{
    RawReadReader* query_reader = RawReadReaderNew(db_dir, INIT_QUERY_DB_TITLE, TRUE);
    CSeqDBInfo volinfo = seqdb_load_volume_info(db_dir, INIT_QUERY_DB_TITLE, vid + 1);
    for (int i = 0; i < volinfo.num_seqs; ++i) {
        RawReadReader_SetReadFlag(query_reader, volinfo.seq_start_id + i);
    }
    RawReadReader_LoadFlaggedReads(query_reader);
    return query_reader;
}

static void
process_one_volume_result(WorkData** wd_array, int num_threads, const char* db_dir, int qi)
{
    if (query_vol_sv_read_is_created(g_wrk_dir, qi)) {
        HBN_LOG("Query volume %d mapping results are processed. Skip it.", qi);
        return;
    }
    char job_name[256];
    sprintf(job_name, "process volume %8d", qi);
    hbn_timing_begin(job_name);
    M4Record* m4a = NULL;
    size_t m4c = 0;
    load_volume_mapping_results(db_dir, qi, 0, &m4a, &m4c);
    HBN_LOG("load %zu m4", m4c);
    RawReadReader* query_reader = s_load_query_reader(db_dir, qi);
    int m4_idx = 0;
    pthread_mutex_t m4_idx_lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_t job_ids[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        wd_array[i]->query_reader = query_reader;
        wd_array[i]->query_vol_index = qi;
        wd_array[i]->m4a = m4a;
        wd_array[i]->m4c = m4c;
        wd_array[i]->m4_idx = &m4_idx;
        wd_array[i]->m4_idx_lock = &m4_idx_lock;
        kv_clear(wd_array[i]->indel_sv_read_list);
        pthread_create(job_ids + i, NULL, thread_func, wd_array[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(job_ids[i], NULL);
        wd_array[i]->query_reader = NULL;
        wd_array[i]->query_vol_index = -1;
        wd_array[i]->m4a = NULL;
        wd_array[i]->m4c = 0;
        wd_array[i]->m4_idx = NULL;
        wd_array[i]->m4_idx_lock = NULL;
    }
    dump_query_vol_sv_read_results(wd_array, num_threads, qi);
    free(m4a);
    RawReadReaderFree(query_reader);
    hbn_timing_end(job_name);
    query_vol_sv_read_make_created(g_wrk_dir, qi);
}

static void
merge_ref_pos_results(int num_query_vols, SvRead* svr_array)
{
    char path[HBN_MAX_PATH_LEN];
    for (int i = 0; i < num_query_vols; ++i) {
        hbn_assert(query_vol_sv_read_is_created(g_wrk_dir, i));
        make_query_vol_sv_read_path(g_wrk_dir, i, path);
        int svrc = 0;
        SvRead* svra = load_sv_read_array(path, &svrc);
        HBN_LOG("merge %d sv reads from volume %d", svrc, i);
        for (int p = 0; p < svrc; ++p) {
            SvRead svr = svra[p];
            svr_array[svr.query_id] = svr;
        }
        free(svra);
    } 
}

static void
print_usage(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s wrk_dir db_dir min_seq_size min_sve_perc_identity max_overhang num_threads [trf_file_path]\n", pn);
}

static void
dump_sv_read_array(WorkData* data, const char* db_dir, SvRead* sv_read_array)
{
    char path[HBN_MAX_PATH_LEN];
    int seqs;
    size_t res;
    char buf1[64], buf2[64];
    CSeqDBInfo dbinfo = seqdb_load_volume_info(db_dir, INIT_QUERY_DB_TITLE, 0);
    const int num_queries = dbinfo.num_seqs;
    CSeqInfo* seqinfo_array = load_seq_infos(db_dir, INIT_QUERY_DB_TITLE, 0, num_queries);
    char* query_names = load_seq_headers(db_dir, INIT_QUERY_DB_TITLE, dbinfo.hdr_offset_from, dbinfo.hdr_offset_to);
    for (int i = 0; i < data->subject_reader->dbinfo.num_seqs; ++i) {
        make_subject_sv_read_path(g_wrk_dir, i, path);
        hbn_dfopen(out, path, "w");
        u64_to_string_datasize(RawReadReader_ReadSize(data->subject_reader, i), buf1);
        HBN_LOG("dumpping queries for subject %d (length: %s)", i, buf1);
        seqs = 0;
        res = 0;
        for (int k = 0; k < num_queries; ++k) {
            SvRead* svr = sv_read_array + k;
            if (svr->subject_id == i) {
                ++seqs;
                res += seqinfo_array[svr->query_id].seq_size;
                const char* name = query_names + seqinfo_array[svr->query_id].hdr_offset;
                dump_sv_read(fprintf, out, *svr, name);
            }
        }
        u64_to_string_comma(seqs, buf1);
        u64_to_string_datasize(res, buf2);
        size_t subject_length = RawReadReader_ReadSize(data->subject_reader, i);
        double cov = 1.0 * res / subject_length;
        HBN_LOG("%s queries (%s, %.2lfx)", buf1, buf2, cov);
        hbn_fclose(out);
    }
    free(query_names);
    free(seqinfo_array);
}

int main(int argc, char* argv[])
{
    if (argc != 7 && argc != 8) {
        print_usage(argv[0]);
        return 1;
    }

    g_wrk_dir = argv[1];
    g_db_dir = argv[2];
    g_min_seq_size = atoi(argv[3]);
    g_min_ident_perc = atof(argv[4]);
    g_max_overhang = atoi(argv[5]);
    const int num_threads = atoi(argv[6]);
    const char* trf_path = (argc == 8) ? argv[7] : NULL;
    create_sv_read_dir(g_wrk_dir);

    pthread_mutex_init(&gq_lock, NULL);
    RawReadReader* subject_reader = RawReadReaderNew(g_db_dir, INIT_SUBJECT_DB_TITLE, FALSE);
    const int num_queries = seqdb_load_num_reads(g_db_dir, INIT_QUERY_DB_TITLE);
    SvRead* sv_read_array = (SvRead*)calloc(num_queries, sizeof(SvRead));
    for (int i = 0; i < num_queries; ++i) init_sv_read(sv_read_array[i]);
    TrfArray* trf_array = trf_path ? TrfArrayBuild(trf_path, g_db_dir, subject_reader->raw_read_info_array) : NULL;

    WorkData* wd_array[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        wd_array[i] = WorkDataNew();
        wd_array[i]->subject_reader = subject_reader;
        wd_array[i]->trf_array = trf_array;
    }
    int num_vols = seqdb_load_num_volumes(g_db_dir, INIT_QUERY_DB_TITLE);
    for (int i = 0; i < num_vols; ++i) {
        process_one_volume_result(wd_array, num_threads, g_db_dir, i);
    }

    merge_ref_pos_results(num_vols, sv_read_array);
    dump_sv_read_array(wd_array[0], g_db_dir, sv_read_array);
    dump_subject_sv_read_file_count(g_wrk_dir, subject_reader->dbinfo.num_seqs);
    for (int i = 0; i < num_threads; ++i) wd_array[i] = WorkDataFree(wd_array[i]);
    free(sv_read_array);
    RawReadReaderFree(subject_reader);
    if (trf_array) trf_array = TrfArrayFree(trf_array);
    return 0;
}
