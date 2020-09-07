#include "chain_and_extend_kmer_matches.h"

#include "hbn_find_subseq_hit.h"
#include "hbn_extend_subseq_hit.h"
#include "hbn_extend_subseq_hit_diff.h"
#include "mecat_results.h"
#include "../../algo/sort_sr_hit_seeds.h"
#include "../../corelib/m4_record.h"
#include "../../corelib/collect_km_mode.h"
#include "../../ncbi_blast/setup/hsp2string.h"

#include <pthread.h>

static const int kExtnQueryBatchSize = 5;

typedef struct {
    CSeqDB* query_vol;
    CSeqDB* subject_vol;
    HbnWordFindData* word_data;
    const HbnProgramOptions* opts;
    int qidx_from;
    int qidx_to;
    int* next_qidx;
    pthread_mutex_t* qidx_lock;
    FILE* out;
    FILE* backup_out;
    pthread_mutex_t* out_lock;
    ChainWorkData* chain;
    HbnSubseqHitExtnData* extn_data;
} ChainAndExtendKmerMatchThreadData;

ChainAndExtendKmerMatchThreadData*
ChainAndExtendKmerMatchThreadDataNew(
    CSeqDB* query_vol,
    CSeqDB* subject_vol,
    HbnWordFindData* word_data,
    const HbnProgramOptions* opts,
    int qidx_from,
    int qidx_to,
    int* next_qidx,
    pthread_mutex_t* qidx_lock,
    FILE* out,
    FILE* backup_out,
    pthread_mutex_t* out_lock)
{
    ChainAndExtendKmerMatchThreadData* data = 
        (ChainAndExtendKmerMatchThreadData*)calloc(1, sizeof(ChainAndExtendKmerMatchThreadData));
    data->query_vol = query_vol;
    data->subject_vol = subject_vol;
    data->word_data = word_data;
    data->opts = opts;
    data->qidx_from = qidx_from;
    data->qidx_to = qidx_to;
    data->next_qidx = next_qidx;
    data->qidx_lock = qidx_lock;
    data->out = out;
    data->backup_out = backup_out;
    data->out_lock = out_lock;
    int min_chain_score = opts->min_ddfs * opts->kmer_size * 0.8;
    data->chain = ChainWorkDataNew(opts->min_ddfs, min_chain_score);
    data->extn_data = HbnSubseqHitExtnDataNew(opts);
    return data;
}

ChainAndExtendKmerMatchThreadData*
ChainAndExtendKmerMatchThreadDataFree(ChainAndExtendKmerMatchThreadData* data)
{
    ChainWorkDataFree(data->chain);
    HbnSubseqHitExtnDataFree(data->extn_data);
    free(data);
    return NULL;
}

static void
fill_init_hit_info(int qid, int qdir, int qsize, int sid, int ssize, HbnInitHit* hit_array, int hit_count)
{
    for (int i = 0; i < hit_count; ++i) {
        HbnInitHit* hit = hit_array + i;
        hit->qid = qid;
        hit->qdir = qdir;
        hit->qsize = qsize;
        hit->sid = sid;
        hit->sdir = FWD;
        hit->ssize = ssize;
        ChainSeed* fcs = hit->chain_seed_array;
        ChainSeed* lcs = fcs + hit->chain_seed_count - 1;
        hit->qbeg = fcs->qoff;
        hit->sbeg = fcs->soff;
        hit->qend = lcs->qoff + lcs->length;
        hit->send = lcs->soff + lcs->length;
        hit->qoff = hit->qend;
        hit->soff = hit->send;
    }
}

static void
find_candidates_for_one_context_one_subject(DDFKmerMatch* ddfkm_array, 
    int ddfkm_count, 
    ChainWorkData* chain, 
    const int kmer_size,
    const int qid,
    const int qdir,
    const int qsize,
    const int sid,
    const idx sbeg,
    const idx send,
    const idx ssize,
    vec_init_hit* hit_list)
{
    kv_clear(chain->fwd_seeds);
    int i = 0;
    while (i < ddfkm_count) {
        int j = i + 1;
        while (j < ddfkm_count && ddfkm_array[i].soff == ddfkm_array[j].soff) ++j;
        if (j - i > 4) {
            i = j;
            continue;
        }
        for (int k = i; k < j; ++k) {
            ChainSeed cs;
            cs.qoff = ddfkm_array[k].qoff;
            cs.sdir = FWD;
            cs.soff = ddfkm_array[k].soff - sbeg;
            cs.length = kmer_size;
            kv_push(ChainSeed, chain->fwd_seeds, cs);            
        }
        i = j;
    }

    ChainSeed* csa = kv_data(chain->fwd_seeds);
    int csc = kv_size(chain->fwd_seeds);
    kv_dinit(vec_init_hit, l_hit_list);
    kv_dinit(vec_chain_seed, l_chain_seed_list);
    chaining_find_candidates(chain, csa, csc, FALSE, FWD, &l_hit_list, &l_chain_seed_list);
    HbnInitHit* hit_array = kv_data(l_hit_list);
    int hit_count = kv_size(l_hit_list);
    fill_init_hit_info(qid, qdir, qsize, sid, ssize, hit_array, hit_count);
    kv_push_v(HbnInitHit, *hit_list, hit_array, hit_count);
    kv_destroy(l_hit_list);
    kv_destroy(l_chain_seed_list);
}

static void
find_candidates_for_one_context(const text_t* db,
    DDFKmerMatch* ddfkm_array, 
    int ddfkm_count, 
    ChainWorkData* chain, 
    const int kmer_size,
    const int qid,
    const int qdir,
    const int qsize,
    vec_init_hit* hit_list)
{
    if (!ddfkm_count) return;

    sort_ddfkm_soff_lt(ddfkm_count, ddfkm_array);

    int i = 0;
    while (i < ddfkm_count) {
        int sid = seqdb_offset_to_seq_id(db, ddfkm_array[i].soff);
        idx sbeg = seqdb_seq_offset(db, sid);
        idx ssize = seqdb_seq_size(db, sid);
        idx send = sbeg + ssize;
        int j = i + 1;
        while (j < ddfkm_count && ddfkm_array[j].soff < send) ++j;
        int n = j - i;
        if (n < chain->min_cnt) {
            i = j;
            continue;
        }
        find_candidates_for_one_context_one_subject(ddfkm_array + i,
            j - i,
            chain,
            kmer_size,
            qid,
            qdir,
            qsize,
            sid,
            sbeg,
            send,
            ssize,
            hit_list);
        i = j;
    }
}

static void
process_one_query(ChainAndExtendKmerMatchThreadData* data,
    BlastHitList* hsp_hit_list,
    vec_subseq_hit* subseq_hit_sink,
    HbnHSPResults* results, int qid)
{
    kv_dinit(vec_init_hit, hit_list);
    CSeqDB* qdb = data->query_vol;
    CSeqDB* sdb = data->subject_vol; 
    hbn_assert(qid < qdb->dbinfo.num_seqs);
    DDFKmerMatch* kma = data->word_data->kmi_array[qid].km_array;
    const int qsize = seqdb_seq_size(qdb, qid);
    int kmc = data->word_data->kmi_array[qid].km_count;
    //HBN_LOG("query %d: kmc = %d", qid, kmc);
    int i = 0;
    while (i < kmc) {
        int j = i + 1;
        while (j < kmc && kma[i].context == kma[j].context) ++j;
        int qdir = kma[i].context&1;
        find_candidates_for_one_context(sdb, kma + i, j - i, data->chain, data->opts->kmer_size, qid, qdir, qsize, &hit_list);
        i = j;
    }
    if (kv_empty(hit_list)) return;

    //HBN_LOG("query %d, hits: %d", qid, kv_size(hit_list));

    HbnSubseqHitExtnData* extn_data = data->extn_data;
    vec_subseq_hit* fwd_sbjct_subseq_list = &extn_data->fwd_sbjct_subseq_list;
    vec_subseq_hit* rev_sbjct_subseq_list = &extn_data->rev_sbjct_subseq_list;
    vec_subseq_hit* sbjct_subseq_list = &extn_data->sbjct_subseq_list;
    const HbnProgramOptions* opts = data->opts;
    kv_dinit(vec_u8, fwd_query_v);
    kv_dinit(vec_u8, rev_query_v);
    HbnInitHit* init_hit_array = kv_data(hit_list);
    int init_hit_count = kv_size(hit_list);
    const char* query_name = seqdb_seq_name(qdb, qid);
    seqdb_extract_sequence(qdb, qid, FWD, &fwd_query_v);
    const u8* fwd_query = kv_data(fwd_query_v);
    seqdb_extract_sequence(qdb, qid, REV, &rev_query_v);
    const u8* rev_query = kv_data(rev_query_v);
    find_candidate_subject_subseqs(init_hit_array,
        init_hit_count,
        qid,
        qsize,
        opts,
        sdb,
        fwd_sbjct_subseq_list,
        rev_sbjct_subseq_list,
        sbjct_subseq_list);  
    if (opts->outfmt >= eSeqid && opts->outfmt <= eSubseqx) {
        kv_push_v(HbnSubseqHit, 
            *subseq_hit_sink, 
            kv_data(*sbjct_subseq_list), 
            kv_size(*sbjct_subseq_list));
    } else {
        hbn_extend_query_subseq_hit_list(kv_data(*sbjct_subseq_list),
                kv_size(*sbjct_subseq_list),
                query_name,
                fwd_query,
                rev_query,
                qid,
                qsize,
                sdb,
                opts,
                extn_data,
                hsp_hit_list,
                results);    
    }

    kv_destroy(hit_list);
    kv_destroy(fwd_query_v);
    kv_destroy(rev_query_v);
}

static void
dump_subseq_hits(HbnSubseqHit* subseq_array,
    const int subseq_count,
    const text_t* query_vol,
    const text_t* subject_vol,
    EOutputFormat outfmt,
    kstring_t* out_buf)
{
    ks_clear(*out_buf);
    for (int i = 0; i < subseq_count; ++i) {
        HbnSubseqHit* subseq = subseq_array + i;
        subseq->qid += query_vol->dbinfo.seq_start_id;
        subseq->sid += subject_vol->dbinfo.seq_start_id;
    }

    char line[1024];
    for (int i = 0; i < subseq_count; ++i) {
        HbnSubseqHit* subseq = subseq_array + i;
        if (outfmt == eSeqid || outfmt == eSeqidx) {
            HbnConsensusInitHit cns_hit;
            cns_hit.qid = subseq->qid;
            cns_hit.sid = subseq->sid;
            cns_hit.score = subseq->score;
            if (outfmt == eSeqid) {
                dump_cns_hit(sprintf, line, cns_hit);
                kputsn(line, strlen(line), out_buf);
            } else {
                const char* input = (const char*)(&cns_hit);
                const int input_len = sizeof(HbnConsensusInitHit);
                kputsn(input, input_len, out_buf);
            }
            continue;
        }
        if (outfmt == eSubseq) {
            dump_subseq_hit(sprintf, line, *subseq);
            kputsn(line, strlen(line), out_buf);
        } else {
            hbn_assert(outfmt == eSubseqx);
            const char* input = (const char*)(subseq);
            const int input_len = sizeof(HbnSubseqHit);
            kputsn(input, input_len, out_buf);
        }
    }
}

static void
dump_m4_hits(const text_t* query_vol,
    const text_t* subject_vol,
    HbnHSPResults* results,
    const HbnProgramOptions* opts,
    const int qidx_from,
    const int qidx_to)
{
    EOutputFormat outfmt = opts->outfmt;
    ks_dinit(line);
    ks_clear(results->output_buf);
    for (int i = qidx_from; i < qidx_to; ++i) {
        BlastHitList* hit_list = results->hitlist_array + i;
        if (hit_list->hsplist_array == 0) continue;
        for (int j = 0; j < hit_list->hsplist_count; ++j) {
            BlastHSPList* hsp_list = hit_list->hsplist_array[j];
            for (int k = 0; k < hsp_list->hspcnt; ++k) {
                BlastHSP* hsp = hsp_list->hsp_array[k];
                hbn_assert(hsp->hbn_subject.strand == FWD);
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
                const char* qname = seqdb_seq_name(query_vol, hsp->hbn_query.oid);
                const char* sname = seqdb_seq_name(subject_vol, hsp->hbn_subject.oid);
                hsp->hbn_query.oid = query_vol->dbinfo.seq_start_id + hsp->hbn_query.oid;
                hsp->hbn_subject.oid = subject_vol->dbinfo.seq_start_id + hsp->hbn_subject.oid;

                if (outfmt >= eM4 && outfmt <= eM4x) {
                    print_one_m4_result(hsp,
                        &results->aligned_strings,
                        qname,
                        sname,
                        opts->dump_cigar,
                        opts->dump_md,
                        outfmt == eM4x,
                        &line,
                        &results->output_buf);
                } else if (outfmt == ePaf) {
                    print_one_paf_result(hsp,
                        &results->aligned_strings,
                        qname,
                        sname,
                        opts->dump_cigar,
                        opts->dump_md,
                        &results->output_buf);
                } else if (outfmt == eSAM) {
                    print_one_sam_result(hsp,
                        &results->aligned_strings,
                        qname,
                        sname,
                        opts->dump_md,
                        &results->output_buf);
                }
            }
        }
    }
    ks_destroy(line);
}

void*
chain_and_extend_km_thread(void* params)
{
    ChainAndExtendKmerMatchThreadData* data = (ChainAndExtendKmerMatchThreadData*)(params);
    const HbnProgramOptions* opts = data->opts;
    HbnHSPResults* results = HbnHSPResultsNew(kExtnQueryBatchSize);
    kv_dinit(vec_subseq_hit, subseq_hit_sink);
    while (1) {
        int qidx_from = 0;
        pthread_mutex_lock(data->qidx_lock);
        qidx_from = *data->next_qidx;
        *data->next_qidx += kExtnQueryBatchSize;
        pthread_mutex_unlock(data->qidx_lock);
        if (qidx_from >= data->qidx_to) break;
        int qidx_to = hbn_min(qidx_from + kExtnQueryBatchSize, data->qidx_to);

        HbnHSPResultsClear(results, 0, kExtnQueryBatchSize);
        kv_clear(subseq_hit_sink);
        for (int i = qidx_from; i < qidx_to; ++i) {
            BlastHitList* hsp_hit_list = results->hitlist_array + (i - qidx_from);
            process_one_query(data, hsp_hit_list, &subseq_hit_sink, results, i);
        }

        if (opts->outfmt >= eSeqid && opts->outfmt <= eSubseqx) {
            dump_subseq_hits(kv_data(subseq_hit_sink),
                kv_size(subseq_hit_sink),
                data->query_vol,
                data->subject_vol,
                opts->outfmt,
                &results->output_buf);
        } else {
            dump_m4_hits(data->query_vol,
                data->subject_vol,
                results,
                opts,
                0,
                kExtnQueryBatchSize);
        }
        pthread_mutex_lock(data->out_lock);
        hbn_fwrite(ks_s(results->output_buf), 1, ks_size(results->output_buf), data->out);
        if (data->backup_out) hbn_fwrite(ks_s(results->output_buf), 1, ks_size(results->output_buf), data->backup_out);
        pthread_mutex_unlock(data->out_lock);
    }
    kv_destroy(subseq_hit_sink);
    HbnHSPResultsFree(results);
    return NULL;
}

void
chain_and_extend_km_mt(HbnWordFindData* word_data,
    const HbnProgramOptions* opts,
    CSeqDB* query_vol,
    CSeqDB* subject_vol,
    const int qidx_from,
    const int qidx_to,
    FILE* out,
    FILE* backup_out)
{
    const int num_threads = opts->num_threads;
    pthread_t jobids[num_threads];
    ChainAndExtendKmerMatchThreadData* data_array[num_threads];
    int next_qidx = qidx_from;
    pthread_mutex_t qidx_lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = ChainAndExtendKmerMatchThreadDataNew(query_vol,
                            subject_vol,
                            word_data,
                            opts,
                            qidx_from,
                            qidx_to,
                            &next_qidx,
                            &qidx_lock,
                            out,
                            backup_out,
                            &out_lock);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(&jobids[i], NULL, chain_and_extend_km_thread, data_array[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobids[i], NULL);
    }
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = ChainAndExtendKmerMatchThreadDataFree(data_array[i]);
    }
}