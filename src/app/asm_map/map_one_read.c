#include "map_one_read.h"

#include "../../algo/hbn_traceback_aux.h"
#include "../../algo/refine_align.h"
#include "../../ncbi_blast/setup/blast_hits.h"
#include "../../ncbi_blast/setup/hsp2string.h"
#include "../map/mecat_results.h"
#include "resovle_broken_map.h"

#define E 100

static BOOL 
is_correct_ovlp(const int qb, const int qe, const int qs, 
    const int sb, const int se, const int ss)
{
    if (qb <= E && qs - qe <= E) return TRUE;
    if (sb <= E && ss - se <= E) return TRUE;
    if (ss - se <= E && qb <= E) return TRUE;
    if (qs - qe <= E && sb <= E) return TRUE;
    return FALSE;
}

static BOOL
ovlp_is_overhang(const int qb, const int qe, const int qs, 
    const int sb, const int se, const int ss)
{
    if (is_correct_ovlp(qb, qe, qs, sb, se, ss)) return FALSE;
    
    BOOL r = (qb <= E) || (qs - qe <= E) || (sb <= E) || (ss - se <= E);
    return r;
}

static void
s_set_blasthsp(const int sr_id,
    const char* sr_name,
    const int sr_dir,
    const int sr_off,
    const int sr_send,
    const int sr_size,
    const int sbjct_id,
    const char* sbjct_name,
    const int sbjct_off,
    const int sbjct_end,
    const int sbjct_size,
    const double perc_identity,
    const int chain_score,
    const int align_score,
    const char* sr_align,
    const char* sbjct_align,
    const int align_size,
    BlastHSP* hsp,
    kstring_t* align_strings)
{
    hsp->hbn_gap_info = NULL;
    hsp->hbn_query.oid = sr_id;
    hsp->hbn_query.offset = sr_off;
    hsp->hbn_query.end = sr_send;
    hsp->hbn_query.strand = sr_dir;
    hsp->hbn_query.seq_size = sr_size;
    hsp->hbn_subject.oid = sbjct_id;
    hsp->hbn_subject.offset = sbjct_off;
    hsp->hbn_subject.end = sbjct_end;
    hsp->hbn_subject.strand = FWD;
    hsp->hbn_subject.seq_size = sbjct_size;
    hsp->hsp_info.perc_identity = perc_identity;
    hsp->hsp_info.raw_score = align_score;
    hsp->hsp_info.chain_score = chain_score;
    hsp->hsp_info.ddf_score = chain_score;
    hsp->hsp_info.query_align_offset = ks_size(*align_strings);
    kputsn(sr_align, align_size, align_strings);
    hsp->hsp_info.subject_align_offset = ks_size(*align_strings);
    kputsn(sbjct_align, align_size, align_strings);
}

static BOOL
chain_seed_list_is_contained(const BlastHSP* hsp_array,
    const int hsp_count,
    const int qdir,
    const int qsize,
    const ChainSeed* chain_seed_array,
    const int chain_seed_count)
{
    const ChainSeed* fcs = chain_seed_array;
    const ChainSeed* lcs = chain_seed_array + chain_seed_count - 1;
    int qbeg = fcs->qoff;
    int qend = lcs->qoff + lcs->length;
    int sbeg = fcs->soff;
    int send = lcs->soff + lcs->length;
    for (int i = 0; i < hsp_count; ++i) {
        const BlastHSP* hsp = hsp_array + i;
        hbn_assert(hsp->hbn_subject.strand == FWD);
        if (hsp->hbn_query.strand != qdir) continue;
        const int X = 200;
        int r = (qbeg + X >= hsp->hbn_query.offset)
                &&
                (qend <= hsp->hbn_query.end + X)
                &&
                (sbeg + X >= hsp->hbn_subject.offset)
                &&
                (send <= hsp->hbn_subject.end + X);
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
    const int X = 200;
    for (int i = 0; i < hsp_count; ++i) {
        BlastHSP* hi = hsp_array + i;
        if (hi->hbn_query.oid == -1) continue;
        for (int j = i + 1; j < hsp_count; ++j) {
            BlastHSP* hj = hsp_array + j;
            if (hj->hbn_query.oid == -1) continue;
            hbn_assert(hj->hsp_info.raw_score <= hi->hsp_info.raw_score);
            if (hj->hbn_query.strand != hi->hbn_query.strand) continue;
            int r = (hj->hbn_query.offset + X >= hi->hbn_query.offset)
                    &&
                    (hj->hbn_query.end <= hi->hbn_query.end + X)
                    &&
                    (hj->hbn_subject.offset + X >= hi->hbn_subject.offset)
                    &&
                    (hj->hbn_subject.end <= hi->hbn_subject.end + X);
            if (r) hj->hbn_query.oid = -1;
        }
    }
    int n = 0;
    for (int i = 0; i < hsp_count; ++i) if (hsp_array[i].hbn_query.oid >= 0) hsp_array[n++] = hsp_array[i];
    *hsp_count_ = n;
}

static int
s_asm_map_one_support_read(const int sr_id,
    const char* sr_name,
    const u8* fwd_sr,
    const u8* rev_sr,
    const int sr_size,
    const int sbjct_id,
    const char* sbjct_name,
    const u8* fwd_sbjct,
    const int sbjct_size,
    const HbnProgramOptions* opts,
    vec_chain_seed* chain_seed_list,
    InitHitFindData* hit_finder,
    HbnTracebackData* tbck_data,
    kstring_t* qaln,
    kstring_t* saln,
    AsmMapResults* results,
    BlastHSPList* hsp_list)
{
    InitHitFindData_AddQuery(hit_finder, sr_id, sr_name, fwd_sr, rev_sr, sr_size);
    InitHitFindData_FindHits(hit_finder, F_R);
    if (kv_empty(hit_finder->hit_list)) return 0;
    const int kMaxExtendHit = opts->max_hsps_per_subject + 2;
    BlastHSP hsp_array[kMaxExtendHit];
    BOOL has_correct_hsp = FALSE;
    BOOL has_overhang_hsp = FALSE;
    int hspcnt = 0;
    for (int i = 0; i < kv_size(hit_finder->hit_list) && i < kMaxExtendHit; ++i) {
        HbnInitHit* hit = &kv_A(hit_finder->hit_list, i);
        //HBN_LOG("extending hit %d / %zu, score = %d", i, kv_size(hit_finder->hit_list), hit->score);
        InitHitFindData_SetupCnsAlignInfoFromHit(hit_finder, hit, chain_seed_list);
        ChainSeed* csa = kv_data(*chain_seed_list);
        const int csc = kv_size(*chain_seed_list);
        if (csc == 0) {
            //HBN_LOG("%d:%s v.s. %d:%s, score = %d, no chain seeds.", sr_id, sr_name, sbjct_id, sbjct_name, hit->score);
            continue;
        }
        if (chain_seed_list_is_contained(hsp_array, hspcnt, hit->sdir, sr_size, csa, csc)) {
            continue;
        }
        const u8* sr = (hit->sdir == FWD) ? fwd_sr : rev_sr;
        validate_mem(HBN_LOG_ARGS_DEFAULT, sr, fwd_sbjct, csa, csc);
        int qb, qe, sb, se;
        double perc_identity;
        int r;
        r = hbn_traceback(tbck_data,
                sr,
                sr_size,
                fwd_sbjct,
                sbjct_size,
                csa,
                csc,
                100,
                0.0,
                TRUE);
        hit->sid = -1;
        if (!r) {
            //HBN_LOG("%d:%s v.s. %d:%s, score = %d, extension fail.", sr_id, sr_name, sbjct_id, sbjct_name, hit->score);
            continue;            
        }
        qb = tbck_data->qoff;
        qe = tbck_data->qend;
        sb =tbck_data->soff;
        se = tbck_data->send;
        perc_identity = tbck_data->ident_perc;
        BlastHSP* hsp = hsp_array + hspcnt;
        s_set_blasthsp(sr_id,
            sr_name,
            hit->sdir,
            qb,
            qe,
            sr_size,
            sbjct_id,
            sbjct_name,
            sb,
            se,
            sbjct_size,
            perc_identity,
            hit->score,
            tbck_data->score,
            tbck_data->qas,
            tbck_data->sas,
            tbck_data->qae - tbck_data->qas,
            hsp,
            &results->align_string);  
        //dump_blasthsp(fprintf, stderr, *hsp);
        if (is_correct_ovlp(qb, qe, sr_size, sb, se, sbjct_size) && perc_identity >= opts->perc_identity) {
            if (hspcnt > 0) {
                //HBN_LOG("correct hsp at %d", i);
                //for (int x = 0; x < hspcnt; ++x) {
                //    dump_blasthsp(fprintf, stderr, hsp_array[x]);
                //}
                //dump_blasthsp(fprintf, stderr, *hsp);
                memcpy(hsp_array, hsp, sizeof(BlastHSP));
            }
            hspcnt = 1;
            has_correct_hsp = TRUE;
            break;
        }
        if (ovlp_is_overhang(qb, qe, sr_size, sb, se, sbjct_size) && perc_identity >= 92.0) {
#if 0
            HBN_LOG("overhang hsp: ");
            dump_blasthsp(fprintf, stderr, *hsp);
            for (size_t x = 0; x < kv_size(hit_finder->hit_list) && x < 6; ++x) {
                hit = &kv_A(hit_finder->hit_list, x);
                print_init_hit_range(hit);
            }
#endif
            has_overhang_hsp =TRUE;
        }
        ++hspcnt;
        if (hspcnt == opts->max_hsps_per_subject) break;
    }
    if (hspcnt == 0) return hspcnt;

#if 1
    if ((!has_correct_hsp) && has_overhang_hsp) {
        int qdir, qb, qe, sb, se, score;
        double perc_identity, eff_perc_identity;
        BOOL r = FALSE;
        r = resolve_broken_map(tbck_data, hsp_array, hspcnt,
                sr_name,
                fwd_sr,
                rev_sr,
                sr_size,
                sbjct_id,
                sbjct_name,
                fwd_sbjct,
                sbjct_size,
                92.0,
                &qdir,
                &qb,
                &qe,
                &sb,
                &se,
                &score,
                &perc_identity,
                &eff_perc_identity,
                qaln,
                saln);
        if (r) {
            BlastHSP* hsp = hsp_array;
            s_set_blasthsp(sr_id,
                sr_name,
                qdir,
                qb,
                qe,
                sr_size,
                sbjct_id,
                sbjct_name,
                sb,
                se,
                sbjct_size,
                eff_perc_identity,
                score,
                score,
                ks_s(*qaln),
                ks_s(*saln),
                ks_size(*qaln),
                hsp,
                &results->align_string);
            hspcnt = 1;
        }
    }
#endif

    ks_introsort_blasthsp_score_gt(hspcnt, hsp_array);
    remove_contained_hsp(hsp_array, &hspcnt);
    hsp_list->oid = sr_id;
    hsp_list->query_index = sbjct_id;
    hsp_list->hspcnt = hspcnt;
    hsp_list->hbn_best_raw_score = hsp_array[0].hsp_info.raw_score;
    for (int i = 0; i < hspcnt; ++i) {
        BlastHSP* hsp = hsp_list->hsp_array[i];
        memcpy(hsp, hsp_array + i, sizeof(BlastHSP));
    }
    return hspcnt;
}

void
asm_map_one_read(AsmMapThreadData* data, int sr_info_idx)
{
    size_t hit_idx_from = data->sr_info_array[sr_info_idx].sr_idx_from;
    size_t hit_idx_to = data->sr_info_array[sr_info_idx].sr_idx_to;
    HbnConsensusInitHit* sr_hit_array = data->hit_array + hit_idx_from;
    const int sr_hit_count = hit_idx_to - hit_idx_from;
    const int subject_id = sr_hit_array[0].sid;
    hbn_assert(subject_id < data->raw_reads->dbinfo.num_seqs);
    const char* subject_name = RawReadReader_ReadName(data->raw_reads, subject_id);
    const int subject_length = RawReadReader_ReadSize(data->raw_reads, subject_id);
    //if (subject_id != 7) return;
    //if (strcmp(subject_name, "29700")) return;
    const u8* fwd_subject = RawReadReader_ExtractRead(data->raw_reads, subject_id, FWD, &data->fwd_subject);
    const u8* rev_subject = RawReadReader_ExtractRead(data->raw_reads, subject_id, REV, &data->rev_subject);
    InitHitFindData* hit_finder = data->hit_finder;
    InitHitFindData_Init(hit_finder, subject_id, subject_name, fwd_subject, rev_subject, subject_length);
    ks_introsort_cns_hit_score_gt(sr_hit_count, sr_hit_array);
    kstring_t* qaln = &data->qaux;
    kstring_t* saln = &data->saux;
    vec_chain_seed* chain_seed_list = &data->chain_seed_list;
    //HBN_LOG("mapping %d:%s:%d:%d", subject_id, subject_name, subject_length, sr_hit_count);
    AsmMapResultsClear(data->results);
    BlastHSPList* hsplist_array = data->results->hsplist_array;
    int hsplist_count = 0;

    for (int i = 0; i < sr_hit_count && i < data->opts->hitlist_size + 5; ++i) {
        const int sr_id = sr_hit_array[i].qid;
        const int sr_size = RawReadReader_ReadSize(data->raw_reads, sr_id);
        const char* sr_name = RawReadReader_ReadName(data->raw_reads, sr_id);
        //if (sr_id == 58325)
        //HBN_LOG("****** mapping support read %s:%d [%d, %d, %d]", sr_name, sr_size, sr_hit_array[i].qid, sr_hit_array[i].sid, sr_hit_array[i].score);
        const u8* fwd_sr = RawReadReader_ExtractRead(data->raw_reads, sr_id, FWD, &data->fwd_read);
        const u8* rev_sr = RawReadReader_ExtractRead(data->raw_reads, sr_id, REV, &data->rev_read);
        if (s_asm_map_one_support_read(sr_id, sr_name, fwd_sr, rev_sr, sr_size,
            subject_id, subject_name, fwd_subject, subject_length,
            data->opts, chain_seed_list, hit_finder, data->tbck_data, 
            qaln, saln,
            data->results, hsplist_array + hsplist_count)) {
            ++hsplist_count;
            if (hsplist_count == data->opts->hitlist_size) break;
        }
    }
    if (!hsplist_count) return;

    ks_clear(data->out_buf);
    kstring_t* line = qaln;
    for (int i = 0; i < hsplist_count; ++i) {
        hbn_assert(hsplist_array[i].oid >= 0);
        const char* sr_name = RawReadReader_ReadName(data->raw_reads, hsplist_array[i].oid);
        for (int j = 0; j < hsplist_array[i].hspcnt; ++j) {
            BlastHSP* hsp = hsplist_array[i].hsp_array[j];
            if (hsp->hbn_query.strand == REV) {
                int offset = hsp->hbn_query.seq_size - hsp->hbn_query.end;
                int end = hsp->hbn_query.seq_size - hsp->hbn_query.offset;
                hsp->hbn_query.offset = offset;
                hsp->hbn_query.end = end;
            }
            if (hsp->hbn_subject.strand == REV) {
                int offset = hsp->hbn_subject.seq_size - hsp->hbn_subject.end;
                int end = hsp->hbn_subject.seq_size - hsp->hbn_subject.offset;
                hsp->hbn_subject.offset = offset;
                hsp->hbn_subject.end = end;
            }
            print_one_m4_result(hsp, &data->results->align_string, sr_name, subject_name,
                FALSE, FALSE, FALSE, line, &data->out_buf);
        }
        hsplist_array[i].oid = -1;
        hsplist_array[i].hspcnt = 0;
        hsplist_array[i].hbn_best_raw_score = 0;
    }

    pthread_mutex_lock(data->out_lock);
    hbn_fwrite(ks_s(data->out_buf), 1, ks_size(data->out_buf), data->out);
    pthread_mutex_unlock(data->out_lock);
    //HBN_LOG("number of overhang: %d", noh);
    //exit(0);
}