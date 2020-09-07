#include "hbn_extend_subseq_hit_diff.h"

#include "../../ncbi_blast/setup/hsp2string.h"

static BOOL
hsp_is_contained(const BlastHSP* hsp_array,
    const int hsp_count,
    BlastHSP* newhsp)
{
    const int E = 200;
    for (int i = 0; i < hsp_count; ++i) {
        const BlastHSP* hsp = hsp_array + i;
        if (hsp->hbn_query.strand != newhsp->hbn_query.strand) continue;
        int r = (newhsp->hbn_query.offset + E >= hsp->hbn_query.offset)
                &&
                (newhsp->hbn_query.end <= hsp->hbn_query.end + E)
                &&
                (newhsp->hbn_subject.offset + E >= hsp->hbn_subject.offset)
                &&
                (newhsp->hbn_subject.end <= hsp->hbn_subject.end + E);
        if (r) return TRUE;
    }
    return FALSE;
}

static BOOL
subseq_is_contained(const BlastHSP* hsp_array,
    const int hsp_count,
    HbnSubseqHit* hit,
    const int qsize,
    BlastHSP* best_hsp)
{
    const int E = 100;
    int qoff = hit->qoff;
    int soff = hit->soff + hit->sfrom;
    for (int i = 0; i < hsp_count; ++i) {
        const BlastHSP* hsp = hsp_array + i;
        if (hsp->hbn_query.strand != hit->qdir) continue;
        int r = (qoff + E >= hsp->hbn_query.offset)
                &&
                (qoff <= hsp->hbn_query.end + E)
                &&
                (soff + E >= hsp->hbn_subject.offset)
                &&
                (soff <= hsp->hbn_subject.end + E);
        if (r) return TRUE;
    }

    int hit_qoff = hit->qoff;
    if (hit->qdir == REV) hit_qoff = qsize - 1 - hit->qoff;
    int hit_soff = hit->soff + hit->sfrom;
    for (int i = -1; i < hsp_count; ++i) {
        const BlastHSP* hsp = (i == -1) ? best_hsp : (hsp_array + i);
        if (best_hsp->hsp_info.raw_score == -1) continue;
        int hsp_qbeg = hsp->hbn_query.offset;
        int hsp_qend = hsp->hbn_query.end;
        if (hsp->hbn_query.strand == REV) {
            hsp_qbeg = qsize - hsp->hbn_query.end;
            hsp_qend = qsize - hsp->hbn_query.offset;
        }
        int hsp_sbeg = hsp->hbn_subject.offset;
        int hsp_send = hsp->hbn_subject.end;
        int r = (hit->qoff >= hsp_qbeg && hit->qoff <= hsp_qend)
                ||
                (hit_soff >= hsp_sbeg && hit_soff <= hsp_send);
        if (!r) continue;
        if (hit->score < hsp->hsp_info.raw_score * 0.4) return TRUE;
    }
    return FALSE;
}

static void
hbn_extend_subject_subseq_hit_list(HbnSubseqHitExtnData* data,
    const text_t* db,
    HbnSubseqHit* hit_array,
    const int hit_count,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_length,
    const HbnProgramOptions* opts,
    BlastHSPList* hsp_list,
    HbnHSPResults* results,
    BlastHSP* best_hsp)
{
    hsp_list->hbn_best_raw_score = 0;
    hsp_list->hspcnt = 0;
    hsp_list->hsp_array = NULL;
    hsp_list->oid = hit_array[0].sid;
    hsp_list->query_index = query_id;

    ks_introsort_subseq_hit_score_gt(hit_count, hit_array);
    BlastHSP hsp_array[opts->max_hsps_per_subject];
    int hsp_count = 0;
    for (int i = 0; i < hit_count && i < opts->max_hsps_per_subject + 1 && hsp_count < opts->max_hsps_per_subject; ++i) {
        HbnSubseqHit* hit = hit_array + i;
        if (subseq_is_contained(hsp_array, hsp_count, hit, query_length, best_hsp)) continue;
        //fprintf(stderr, "extending \t");
        //dump_subseq_hit(fprintf, stderr, *hit);
        const u8* subject = db->unpacked_seq + seqdb_seq_offset(db, hit->sid) + hit->sfrom;
        const int subject_length = hit->sto - hit->sfrom;
        const u8* query = (hit->qdir == FWD) ? fwd_query : rev_query;
        data->diff_data->qid = query_id;
        data->diff_data->sid = hit->sid;
        int r = diff_align(data->diff_data,
                    query,
                    hit->qoff,
                    query_length,
                    subject,
                    hit->soff,
                    subject_length,
                    opts->query_cov_hsp_res,
                    opts->perc_identity,
                    !opts->skip_overhang);
        if (!r) continue;
        int qoff = data->diff_data->qoff;
        int qend = data->diff_data->qend;
        int soff = data->diff_data->soff;
        int send = data->diff_data->send;
        double ident_perc = data->diff_data->ident_perc;

            BlastHSP* hsp = hsp_array + hsp_count;
            hsp->hbn_gap_info = NULL;
            hsp->hbn_query.oid = query_id;
            hsp->hbn_query.offset = qoff;
            hsp->hbn_query.end = qend;
            hsp->hbn_query.strand = hit->qdir;
            hsp->hbn_query.seq_size = query_length;
            hsp->hbn_subject.oid = hit->sid;
            hsp->hbn_subject.offset = soff + hit->sfrom;
            hsp->hbn_subject.end = send + hit->sfrom;
            hsp->hbn_subject.strand = FWD;
            hsp->hbn_subject.seq_size = seqdb_seq_size(db, hit->sid);
            if (hsp_is_contained(hsp_array, hsp_count, hsp)) continue;
            hsp->hsp_info.perc_identity = ident_perc;
            hsp->hsp_info.raw_score = data->diff_data->score;
            hsp->hsp_info.chain_score = hit->score;
            const char* qaln = data->diff_data->qas;
            const char* saln = data->diff_data->sas;
            const int aln_len = data->diff_data->qae - data->diff_data->qas;
            hsp->hsp_info.query_align_offset = ks_size(results->aligned_strings);
            kputsn(qaln, aln_len, &results->aligned_strings);
            hsp->hsp_info.subject_align_offset = ks_size(results->aligned_strings);
            kputsn(saln, aln_len, &results->aligned_strings);
            ++hsp_count;
            if (hsp->hsp_info.raw_score > best_hsp->hsp_info.raw_score) memcpy(best_hsp, hsp, sizeof(BlastHSP));
            //dump_blasthsp(fprintf, stderr, *hsp);
            if (hsp_count == opts->max_hsps_per_subject) break;
    }
    if (!hsp_count) return;

    ks_introsort_blasthsp_score_gt(hsp_count, hsp_array);
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
hbn_extend_query_subseq_hit_list_diff(HbnSubseqHit* subseq_hit_array,
    int subseq_hit_count,
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
    BlastHSPList hsplist_array[opts->hitlist_size];
    int hsplist_count = 0;
    int i = 0;
    BlastHSP best_hsp; best_hsp.hsp_info.raw_score = -1;
    while (i < subseq_hit_count) {
        int j = i + 1;
        while (j < subseq_hit_count && subseq_hit_array[i].sid == subseq_hit_array[j].sid) ++j;
        BlastHSPList* hsp_list = hsplist_array + hsplist_count;
        hbn_extend_subject_subseq_hit_list(data,
            db,
            subseq_hit_array + i,
            j - i,
            query_id,
            fwd_query,
            rev_query,
            query_length,
            opts,
            hsp_list,
            results,
            &best_hsp);
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