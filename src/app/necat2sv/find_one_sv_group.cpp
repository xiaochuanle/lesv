#include "find_one_sv_group.h"

#include <algorithm>
#include <vector>

using namespace std;

KHASH_SET_INIT_INT(int_set);

static const int kSvsigWindow = 10;
static const int kMinSvsigCnt = 4;
static const int kMaxSvsigDist = 20;

static const int kSvsigWindowRelax = 50;
static const int kMinSvsigCntRelax = 4;
static const int kMaxSvsigDistRelax = 100;
static const int kMaxSvsigLenDiff = 50;
static const double kMaxSvsigLenDiffRatio = 0.1;

static void
s_dump_svsig_group(SvSignature* sga, int sgc, FILE* out, int* group_id, void* added_id_set_)
{
    khash_t(int_set)* added_id_set = (khash_t(int_set)*)(added_id_set_);
    fprintf(out, "%d\t%d\n", *group_id, sgc);
    ++(*group_id);
    int r;
    for (int i = 0; i < sgc; ++i) {
        dump_svsig(fprintf, out, sga[i], NULL);
        sga[i].type = eGapAlignInvalid;
        kh_put(int_set, added_id_set, sga[i].qid, &r);
    }
}

/////////////////////// ins relax

static BOOL
s_on_same_ins_group_relax(SvSignature* lhs, SvSignature* rhs)
{
    int lhs_l = lhs->qto - lhs->qfrom;
    int rhs_l = rhs->qto - rhs->qfrom;
    int max_l = hbn_max(lhs_l, rhs_l);
    int min_l = hbn_min(lhs_l, rhs_l);
    hbn_assert(min_l > 0);
    return ((max_l - min_l) <= (max_l * kMaxSvsigLenDiffRatio))
           &&
           (max_l - min_l <= kMaxSvsigLenDiff);
}

extern "C"
BOOL
find_next_ins_group_relax(SvSignature* sga,
    int sgc,
    int* next_i,
    int* group_id,
    FILE* out,
    void* added_id_set)
{
    int i = *next_i;
    int last_i = *next_i;
    int last_send = sga[last_i].sfrom + kSvsigWindowRelax;
    int cnt = 0;
    int j = last_i + 1;
    kv_dinit(vec_sv_sig, svsig_list);
    while (j < sgc) {
        if (sga[j].sfrom > last_send) break;
        if (s_on_same_ins_group_relax(&sga[last_i], &sga[j])) {
            last_i = j;
            last_send = sga[j].sfrom + kSvsigWindowRelax;
            ++cnt;
            kv_push(SvSignature, svsig_list, sga[j]);
        }
        ++j;
    }

    if (cnt < kMinSvsigCntRelax) {
        kv_destroy(svsig_list);
        *next_i = i + 1;
        return FALSE;
    }

    kv_push(SvSignature, svsig_list, sga[i]);
    s_dump_svsig_group(kv_data(svsig_list), kv_size(svsig_list), out, group_id, added_id_set);

#if 0
    HBN_LOG("find ins relax group %d", *group_id);
    for (int k = 0; k < kv_size(svsig_list); ++k) {
        SvSignature sig = kv_A(svsig_list, k);
        fprintf(stderr, "%d\t", sig.qto - sig.qfrom);
        dump_svsig(fprintf, stderr, sig, NULL);
    }
#endif
    kv_destroy(svsig_list);
    *next_i = last_i + 1;
    return TRUE;    
}

////////////////// end ins relax

extern "C"
BOOL 
find_next_ins_group(SvSignature* sga, 
    int sgc, 
    int* next_i, 
    int* group_id, 
    FILE* out,
    void* added_id_set)
{
    int i = *next_i;
    int j = i + 1;
    int soff = sga[i].sfrom;
    int send = soff + kSvsigWindow;
    int cnt = 0;
    while (j < sgc && sga[j].sfrom < send) ++j;
    cnt = j - i;
    
    int max_cnt = cnt;
    int max_i = i + cnt / 2;
    while (1) {
        if (j >= sgc) break;

        int reduced = 0;
        for (int k = i; k < j; ++k) if (sga[k].sfrom == soff) ++reduced;
        int k = j;
        int added = 0;
        while (k < sgc && sga[k].sfrom == send + 1) {
            ++added;
            ++k;
        }
        if (added == 0) break;

        cnt -= reduced;
        cnt += added;
        ++soff;
        ++send;
        i += reduced;
        j = k;
        if (cnt > max_cnt) {
            max_cnt = cnt;
            max_i = i + cnt / 2;
        }
    }

    int gi_from = max_i;
    while (gi_from > (*next_i) && sga[max_i].sfrom - sga[gi_from-1].sfrom <= kMaxSvsigDist) --gi_from;
    int gi_to = max_i + 1;
    while (gi_to < sgc && sga[gi_to].sfrom - sga[max_i].sfrom <= kMaxSvsigDist) ++gi_to;
    cnt = gi_to - gi_from;
    if (cnt >= kMinSvsigCnt) s_dump_svsig_group(sga + gi_from, gi_to - gi_from, out, group_id, added_id_set);

#if 0
    HBN_LOG("find ins group %d max_i = %d, [%d, %d, %d]", *group_id, max_i, gi_from, gi_to, sgc);
    ++group_id;
    for (int k = gi_from; k < gi_to; ++k) {
        fprintf(stderr, "%d\t%d\t", k, sga[k].qto - sga[k].qfrom);
        dump_svsig(fprintf, stderr, sga[k], NULL);
    }
#endif

    *next_i = gi_to;
    return TRUE;
}

/////////////////////// del relax

static BOOL
s_on_same_del_group_relax(SvSignature* lhs, SvSignature* rhs)
{
    int lhs_l = lhs->sto - lhs->sfrom;
    int rhs_l = rhs->sto - rhs->sfrom;
    int max_l = hbn_max(lhs_l, rhs_l);
    int min_l = hbn_min(lhs_l, rhs_l);
    hbn_assert(min_l > 0);
    return ((max_l - min_l) <= (max_l * kMaxSvsigLenDiffRatio))
           &&
           (max_l - min_l <= kMaxSvsigLenDiff);
}

extern "C"
BOOL
find_next_del_group_relax(SvSignature* sga,
    int sgc,
    int* next_i,
    int* group_id,
    FILE* out,
    void* added_id_set)
{
    int i = *next_i;
    int last_i = *next_i;
    int last_send = sga[last_i].sfrom + kSvsigWindowRelax;
    int cnt = 0;
    int j = last_i + 1;
    kv_dinit(vec_sv_sig, svsig_list);
    while (j < sgc) {
        if (sga[j].sfrom > last_send) break;
        if (s_on_same_del_group_relax(&sga[last_i], &sga[j])) {
            last_i = j;
            last_send = sga[j].sfrom + kSvsigWindowRelax;
            ++cnt;
            kv_push(SvSignature, svsig_list, sga[j]);
        }
        ++j;
    }

    if (cnt < kMinSvsigCntRelax) {
        *next_i = i + 1;
        kv_destroy(svsig_list);
        return FALSE;
    }

    kv_push(SvSignature, svsig_list, sga[i]);
    s_dump_svsig_group(kv_data(svsig_list), kv_size(svsig_list), out, group_id, added_id_set);

#if 0
    HBN_LOG("find del relax group %d", *group_id);
    for (int k = 0; k < kv_size(svsig_list); ++k) {
        SvSignature sig = kv_A(svsig_list, k);
        fprintf(stderr, "%d\t", sig.sto - sig.sfrom);
        dump_svsig(fprintf, stderr, sig, NULL);
    }
#endif

    kv_destroy(svsig_list);
    *next_i = last_i + 1;
    return TRUE;    
}

////////////////// end del relax

extern "C"
BOOL 
find_next_del_group(SvSignature* sga, 
    int sgc, 
    int* next_i, 
    int* group_id, 
    FILE* out,
    void* added_id_set)
{
    int i = *next_i;
    int j = i + 1;
    int soff = sga[i].sfrom;
    int send = soff + kSvsigWindow;
    int cnt = 0;
    while (j < sgc && sga[j].sfrom < send) ++j;
    cnt = j - i;
    
    int max_cnt = cnt;
    int max_i = i + cnt / 2;
    while (1) {
        if (j >= sgc) break;

        int reduced = 0;
        for (int k = i; k < j; ++k) if (sga[k].sfrom == soff) ++reduced;
        int k = j;
        int added = 0;
        while (k < sgc && sga[k].sfrom == send + 1) {
            ++added;
            ++k;
        }
        if (added == 0) break;

        cnt -= reduced;
        cnt += added;
        ++soff;
        ++send;
        i += reduced;
        j = k;
        if (cnt > max_cnt) {
            max_cnt = cnt;
            max_i = i + cnt / 2;
        }
    }

    int gi_from = max_i;
    while (gi_from > (*next_i) && sga[max_i].sfrom - sga[gi_from-1].sfrom <= kMaxSvsigDist) --gi_from;
    int gi_to = max_i + 1;
    while (gi_to < sgc && sga[gi_to].sfrom - sga[max_i].sfrom <= kMaxSvsigDist) ++gi_to;

    cnt = gi_to - gi_from;
    if (cnt >= kMinSvsigCnt) s_dump_svsig_group(sga + gi_from, cnt, out, group_id, added_id_set);

#if 0
    HBN_LOG("find del group %d max_i = %d, [%d, %d, %d]", group_id, max_i, gi_from, gi_to, sgc);
    ++group_id;
    for (int k = gi_from; k < gi_to; ++k) {
        fprintf(stderr, "%d\t%d\t", k, sga[k].sto - sga[k].sfrom);
        dump_svsig(fprintf, stderr, sga[k], NULL);
    }
#endif

    *next_i = gi_to;
    return TRUE;
}