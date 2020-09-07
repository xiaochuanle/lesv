#include "hbn_extend_subseq_hit.h"

#include "../../algo/hbn_traceback_aux.h"
#include "../../ncbi_blast/setup/hsp2string.h"
#include "../../algo/refine_align.h"
#include "../../corelib/ksort.h"

HbnSubseqHitExtnData*
HbnSubseqHitExtnDataNew(const HbnProgramOptions* opts)
{
    HbnSubseqHitExtnData* data = (HbnSubseqHitExtnData*)calloc(1, sizeof(HbnSubseqHitExtnData));
    if (opts->skip_memsc) {
        data->diff_data = DiffGapAlignDataNew();
    } else {
        const int memsc_kmer_size = opts->memsc_kmer_size;
        const int memsc_mem_size = opts->memsc_mem_size;
        hbn_assert(memsc_kmer_size <= memsc_mem_size);
        const int memsc_window_size = opts->memsc_kmer_window;
        const int memsc_score = opts->memsc_score;
        data->hit_finder = InitHitFindDataNew(memsc_kmer_size,
                                memsc_window_size,
                                memsc_score);
        data->traceback_data = HbnTracebackDataNew();
    }
    kv_init(data->chain_seed_list);
    kv_init(data->tbck_chain_seed_list);
    kv_init(data->fwd_sbjct_subseq_list);
    kv_init(data->rev_sbjct_subseq_list);
    kv_init(data->sbjct_subseq_list);
    return data;
}

HbnSubseqHitExtnData*
HbnSubseqHitExtnDataFree(HbnSubseqHitExtnData* data)
{
    if (data->hit_finder) data->hit_finder = InitHitFindDataFree(data->hit_finder);
    if (data->traceback_data) HbnTracebackDataFree(data->traceback_data);
    if (data->diff_data) DiffGapAlignDataFree(data->diff_data);
    kv_destroy(data->chain_seed_list);
    kv_destroy(data->tbck_chain_seed_list);
    kv_destroy(data->fwd_sbjct_subseq_list);
    kv_destroy(data->rev_sbjct_subseq_list);
    kv_destroy(data->sbjct_subseq_list);
    free(data);
    return NULL;
}

static BOOL
subseq_is_contained(HbnSubseqHit* hit_array, int hit_idx)
{
    const int E = 200;
    int sfrom = hit_array[hit_idx].sfrom;
    int sto = hit_array[hit_idx].sto;
    int qdir = hit_array[hit_idx].qdir;
    for (int i = 0; i < hit_idx - 1; ++i) {
        if (hit_array[i].qdir != qdir) continue;
        int ss = hit_array[i].sfrom;
        int st = hit_array[i].sto;
        if (sfrom + E >= ss && sto <= st + E) return TRUE;
    }
    return FALSE;
}

static BOOL
chain_seed_list_is_contained(const BlastHSP* hsp_array,
    const int hsp_count,
    const int qdir,
    const int qsize,
    const int sid,
    const int sfrom,
    const int chain_score,
    const ChainSeed* chain_seed_array,
    const int chain_seed_count)
{
    const ChainSeed* fcs = chain_seed_array;
    const ChainSeed* lcs = chain_seed_array + chain_seed_count - 1;
    int qbeg = fcs->qoff;
    int qend = lcs->qoff + lcs->length;
    int sbeg = fcs->soff + sfrom;
    int send = lcs->soff + lcs->length + sfrom;
    for (int i = 0; i < hsp_count; ++i) {
        const BlastHSP* hsp = hsp_array + i;
        hbn_assert(hsp->hbn_subject.strand == FWD);
        if (hsp->hbn_query.strand != qdir) continue;
        const int E = 200;
        int r = (qbeg + E >= hsp->hbn_query.offset)
                &&
                (qend <= hsp->hbn_query.end + E)
                &&
                (sbeg + E >= hsp->hbn_subject.offset)
                &&
                (send <= hsp->hbn_subject.end + E);
        if (r) {
            return TRUE;
        }
    }

    return FALSE;
}

static void
remove_contained_hsp(BlastHSP* hsp_array, int* hsp_count_)
{
    int hsp_count = *hsp_count_;
    ks_introsort_blasthsp_score_gt(hsp_count, hsp_array);
    const int E = 200;
    for (int i = 0; i < hsp_count; ++i) {
        BlastHSP* hi = hsp_array + i;
        if (hi->hbn_query.oid == -1) continue;
        for (int j = i + 1; j < hsp_count; ++j) {
            BlastHSP* hj = hsp_array + j;
            if (hj->hbn_query.oid == -1) continue;
            hbn_assert(hj->hsp_info.raw_score <= hi->hsp_info.raw_score);
            if (hj->hbn_query.strand != hi->hbn_query.strand) continue;
            int r = (hj->hbn_query.offset + E >= hi->hbn_query.offset)
                    &&
                    (hj->hbn_query.end <= hi->hbn_query.end + E)
                    &&
                    (hj->hbn_subject.offset + E >= hi->hbn_subject.offset)
                    &&
                    (hj->hbn_subject.end <= hi->hbn_subject.end + E);
            if (r) hj->hbn_query.oid = -1;
        }
    }
    int n = 0;
    for (int i = 0; i < hsp_count; ++i) if (hsp_array[i].hbn_query.oid >= 0) hsp_array[n++] = hsp_array[i];
    *hsp_count_ = n;
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
                HBN_LOG("**** find DEL pos at [%d] x [%d, %d), len = %d", qi, si, si + n, n);
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
                HBN_LOG("**** find INS pos at [%d, %d] x [%d], length = %d", qi, qi + n, si, n);
            }
            qi += n;
            i = j;
            continue;
        }
    }
    hbn_assert(qi == qe);
    hbn_assert(si == subject_offset + se);
}

static BOOL
s_is_target_query(const char* query_name)
{
    const char* target1 = "80f02c09-1ecd-4568-a256-93bffd841e80";
    const char* target2 = "20e4c353-6b85-4e8e-86e3-63984260199d";
    const char* pos = strstr(query_name, target1);
    if (pos != NULL) return TRUE;
    pos = strstr(query_name, target2);
    if (pos != NULL) return TRUE;
    return FALSE;
}

static void
s_hbn_extend_subject_subseq_hit_list(HbnSubseqHitExtnData* data,
    const text_t* db,
    HbnSubseqHit* hit_array,
    const int hit_count,
    const int query_id,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_length,
    const HbnProgramOptions* opts,
    BlastHSPList* hsp_list,
    HbnHSPResults* results)
{
    hsp_list->hbn_best_raw_score = 0;
    hsp_list->hspcnt = 0;
    hsp_list->hsp_array = NULL;
    hsp_list->oid = hit_array[0].sid;
    hsp_list->query_index = query_id;

    ks_introsort_subseq_hit_score_gt(hit_count, hit_array);
    InitHitFindData* hit_finder = data->hit_finder;
    BlastHSP hsp_array[opts->max_hsps_per_subject];
    kv_dinit(vec_u8, subject_v);
    int hsp_count = 0;
    //HBN_LOG("%s, hit_count: %d", query_name, hit_count);
    for (int i = 0; i < hit_count && i < opts->max_hsps_per_subject + 2 && hsp_count < opts->max_hsps_per_subject; ++i) {
        HbnSubseqHit* hit = hit_array + i;
        if (subseq_is_contained(hit_array, i)) {
            HBN_LOG("subseq is contained");
            continue;
        }
        //fprintf(stderr, "extending %d\t", i);
        //dump_subseq_hit(fprintf, stderr, *hit);
        seqdb_extract_subsequence(db, hit->sid, hit->sfrom, hit->sto, FWD, &subject_v);
        //const u8* subject = db->unpacked_seq + seqdb_seq_offset(db, hit->sid) + hit->sfrom;
        const u8* subject = kv_data(subject_v);
        const int subject_length = hit->sto - hit->sfrom;
        int strand = (opts->align_task == eHbnTask_pm) ? F_R : hit->qdir;
        InitHitFindData_AddQuery(hit_finder, hit->sid, NULL, subject, NULL, subject_length);
        InitHitFindData_FindHits(hit_finder, strand);
        const int max_k = (opts->align_task == eHbnTask_pm) ? opts->max_hsps_per_subject : 3;
        for (size_t k = 0; k < kv_size(hit_finder->hit_list) && k < max_k; ++k) {
            HbnInitHit* init_hit = kv_data(hit_finder->hit_list) + k;
            InitHitFindData_SetupMapAlignInfoFromHit(hit_finder,
                init_hit,
                &data->chain_seed_list);
            ChainSeed* csa = kv_data(data->chain_seed_list);
            int csc = kv_size(data->chain_seed_list);
            if (csc == 0) continue;

            if (chain_seed_list_is_contained(hsp_array, 
                    hsp_count, 
                    init_hit->sdir,
                    query_length,
                    hit->sid,
                    hit->sfrom,
                    init_hit->score,
                    csa,
                    csc)) {
                //HBN_LOG("chain seed is contained");
                continue;
            }
            //HBN_LOG("qdir = %d, score = %d", init_hit->sdir, init_hit->score);
            if (!hbn_traceback(data->traceback_data,
                    (init_hit->sdir == FWD) ? fwd_query : rev_query,
                    query_length,
                    subject,
                    subject_length,
                    csa,
                    csc,
                    opts->query_cov_hsp_res,
                    opts->perc_identity,
                    !opts->skip_overhang)) {
                //HbnTracebackDataDump(fprintf, stderr, data->traceback_data);
                continue;
            }
            BlastHSP* hsp = hsp_array + hsp_count;
            hsp->hbn_gap_info = NULL;
            hsp->hbn_query.oid = query_id;
            hsp->hbn_query.offset = data->traceback_data->qoff;
            hsp->hbn_query.end = data->traceback_data->qend;
            hsp->hbn_query.strand = hit->qdir;
            hsp->hbn_query.seq_size = query_length;
            hsp->hbn_subject.oid = hit->sid;
            hsp->hbn_subject.offset = data->traceback_data->soff + hit->sfrom;
            hsp->hbn_subject.end = data->traceback_data->send + hit->sfrom;
            hsp->hbn_subject.strand = FWD;
            hsp->hbn_subject.seq_size = seqdb_seq_size(db, hit->sid);
            hsp->hsp_info.perc_identity = data->traceback_data->ident_perc;
            hsp->hsp_info.raw_score = data->traceback_data->score;
            hsp->hsp_info.chain_score = init_hit->score;
            hsp->hsp_info.ddf_score = hit->score;
            const char* qaln = data->traceback_data->qas;
            const char* saln = data->traceback_data->sas;
            const int aln_len = data->traceback_data->qae - data->traceback_data->qas;
            if (0 && s_is_target_query(query_name)) {
                //dump_align_string(qaln, saln, aln_len, stderr);
                ks_dinit(qas);
                ks_dinit(sas);
                int qb, qe, sb, se;
                double perc_identity;
                int r = nw_ksw2_extd2(data->traceback_data->ksw,
                            query_id,
                            (init_hit->sdir == FWD) ? fwd_query : rev_query,
                            data->traceback_data->qoff,
                            data->traceback_data->qend,
                            query_length,
                            hit->sid,
                            subject,
                            data->traceback_data->soff,
                            data->traceback_data->send,
                            subject_length,
                            0,
                            0,
                            hbn_max(query_length, subject_length),
                            &qb,
                            &qe,
                            &sb,
                            &se,
                            &perc_identity,
                            &qas,
                            &sas);
                if (r) {
                    //dump_align_string(ks_s(qas), ks_s(sas), ks_size(qas), stderr);
                    HBN_LOG("[%d, %d, %d] x [%d, %d, %d], %g", qb, qe, query_length,
                        sb, se, subject_length, perc_identity);
                    s_find_sv_signature(query_name, 
                        ks_s(qas),
                        ks_s(sas),
                        ks_size(qas),
                        0,
                        hsp->hbn_query.offset,
                        hsp->hbn_query.end,
                        hsp->hbn_subject.offset,
                        hsp->hbn_subject.end);
                }
                ks_destroy(qas);
                ks_destroy(sas);
            }
            hsp->hsp_info.query_align_offset = ks_size(results->aligned_strings);
            kputsn(qaln, aln_len, &results->aligned_strings);
            hsp->hsp_info.subject_align_offset = ks_size(results->aligned_strings);
            kputsn(saln, aln_len, &results->aligned_strings);
            ++hsp_count;
            //dump_blasthsp(fprintf, stderr, *hsp);
            if (hsp_count == opts->max_hsps_per_subject) break;
        }
    }
    kv_destroy(subject_v);
    if (!hsp_count) return;

    remove_contained_hsp(hsp_array, &hsp_count);
    hsp_list->hbn_best_raw_score = hsp_array[0].hsp_info.raw_score;
    hsp_list->hspcnt = hsp_count;
    hsp_list->hsp_array = (BlastHSP**)SmallObjectAllocAlloc(results->pointer_alloc, hsp_count);
    for (int i = 0; i < hsp_count; ++i) {
        BlastHSP* hsp = (BlastHSP*)SmallObjectAllocAlloc(results->hsp_alloc, 1);
        memcpy(hsp, hsp_array + i, sizeof(BlastHSP));
        hsp_list->hsp_array[i] = hsp;
    }
}

void
hbn_extend_query_subseq_hit_list(HbnSubseqHit* subseq_hit_array,
    int subseq_hit_count,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_id,
    const int query_length,
    const text_t* db,
    const HbnProgramOptions* opts,
    HbnSubseqHitExtnData* data,
    BlastHitList* hit_list,
    HbnHSPResults* results)
{
    InitHitFindData_Init(data->hit_finder,
        query_id,
        query_name,
        fwd_query,
        rev_query,
        query_length);

    BlastHSPList hsplist_array[opts->hitlist_size];
    int hsplist_count = 0;
    int i = 0;
    while (i < subseq_hit_count) {
        int j = i + 1;
        while (j < subseq_hit_count && subseq_hit_array[i].sid == subseq_hit_array[j].sid) ++j;
        BlastHSPList* hsp_list = hsplist_array + hsplist_count;
        s_hbn_extend_subject_subseq_hit_list(data,
            db,
            subseq_hit_array + i,
            j - i,
            query_id,
            query_name,
            fwd_query,
            rev_query,
            query_length,
            opts,
            hsp_list,
            results);
        if (hsp_list->hspcnt) ++hsplist_count;
        i = j;
    }

    if (!hsplist_count) return;
    hit_list->hsplist_array = (BlastHSPList**)SmallObjectAllocAlloc(results->pointer_alloc, hsplist_count);
    hit_list->hsplist_count = hsplist_count;
    hit_list->hsplist_max = hsplist_count;
    for (i = 0; i < hsplist_count; ++i) {
        BlastHSPList* hsp_list = (BlastHSPList*)SmallObjectAllocAlloc(results->hsplist_alloc, 1);
        memcpy(hsp_list, hsplist_array + i, sizeof(BlastHSPList));
        hit_list->hsplist_array[i] = hsp_list;
    }
}