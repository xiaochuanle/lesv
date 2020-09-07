#include "hbn_word_finder.h"

#include "sort_sr_hit_seeds.h"
#include "../corelib/ksort.h"

#include <pthread.h>

static const int kSeedingSeqSize = 300;
static const int kSeedingSeqStride = 200;

WordFinderThreadData*
WordFinderThreadDataNew(
    LookupTable* lktbl,
    const text_t* subject_vol,
    const text_t* query_vol,
    int qidx_from,
    int qidx_to,
    int thread_id,
    int num_threads,
    QueryKmerMatchInfo* kmi_array,
    int kmer_size,
    int kmer_window,
    int max_kmer_occ)
{
    WordFinderThreadData* data = (WordFinderThreadData*)calloc(1, sizeof(WordFinderThreadData));
    data->lktbl = lktbl;
    data->subject_vol = subject_vol;
    data->query_vol = query_vol;
    data->qidx_from = qidx_from;
    data->qidx_to = qidx_to;
    data->thread_id = thread_id;
    data->num_threads = num_threads;
    kv_init(data->l_ddfkm_list);
    data->g_kmi_array = kmi_array;
    data->kmer_size = kmer_size;
    data->kmer_window = kmer_window;
    data->max_kmer_occ = max_kmer_occ;
    return data;
}

WordFinderThreadData*
WordFinderThreadDataFree(WordFinderThreadData* data)
{
    kv_destroy(data->l_ddfkm_list);
    free(data);
    return NULL;
}

HbnWordFindData*
HbnWordFindDataNew(LookupTable* lktbl,
    const text_t* subject_vol,
    const text_t* query_vol,
    int num_threads,
    int kmer_size,
    int kmer_window,
    int max_kmer_occ)
{
    HbnWordFindData* data = (HbnWordFindData*)calloc(1, sizeof(HbnWordFindData));
    data->lktbl = lktbl;
    data->subject_vol = subject_vol;
    data->query_vol = query_vol;
    data->num_threads = num_threads;
    data->wftd_array = (WordFinderThreadData**)calloc(num_threads, sizeof(WordFinderThreadData*));
    data->kmi_array = (QueryKmerMatchInfo*)calloc(query_vol->dbinfo.num_seqs, sizeof(QueryKmerMatchInfo));
    for (int i = 0; i < num_threads; ++i) {
        data->wftd_array[i] = WordFinderThreadDataNew(lktbl,
                                subject_vol,
                                query_vol,
                                0,
                                0,
                                i,
                                num_threads,
                                data->kmi_array,
                                kmer_size,
                                kmer_window,
                                max_kmer_occ);
    }
    data->kmer_size = kmer_size;
    data->kmer_window = kmer_window;
    data->max_kmer_occ = max_kmer_occ;
    return data;
}

HbnWordFindData*
HbnWordFindDataFree(HbnWordFindData* data)
{
    for (int i = 0; i < data->num_threads; ++i) {
        data->wftd_array[i] = WordFinderThreadDataFree(data->wftd_array[i]);
    }
    free(data->wftd_array);
    free(data->kmi_array);
    free(data);
    return NULL;
}

void
HbnWordFindData_FixThreadQidxRange(HbnWordFindData* data, const int qidx_from, const int qidx_to)
{
    int num_threads = data->num_threads;
    IntPair qidx_range[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        qidx_range[i].first = 1;
        qidx_range[i].second = 0;
    }

    const text_t* query_vol = data->query_vol;
    size_t res = 0;
    for (int i = qidx_from; i < qidx_to; ++i) res += seqdb_seq_size(query_vol, i);
    size_t res_pt = (res + num_threads - 1) / num_threads;
    int qid = qidx_from;
    int tid = 0;
    char res_buf[64];
    while (tid < num_threads && qid < qidx_to) {
        size_t curr = seqdb_seq_size(query_vol, qid);
        int next = qid + 1;
        while (next < qidx_to && curr < res_pt) {
            curr += seqdb_seq_size(query_vol, next);
            ++next;
        }
        u64_to_string_datasize(curr, res_buf);
        HBN_LOG("thread %d: query %d --- %d (%d, %s)", tid, qid, next, next - qid, res_buf);
        qidx_range[tid].first = qid;
        qidx_range[tid].second = next;
        ++tid;
        qid = next;
    }
    hbn_assert(num_threads > 0);
    if (tid >= num_threads) tid = num_threads - 1;
    qidx_range[tid].second = qidx_to;

    for (int i = 0; i < num_threads; ++i) {
        data->wftd_array[i]->qidx_from = qidx_range[i].first;
        data->wftd_array[i]->qidx_to = qidx_range[i].second;
    }
}

static int
extract_hash_values(const u8* read,
    const int read_size,
    const int kmer_size,
    const int window_size,
    u64* hash_array)
{
    if (read_size < kmer_size) return 0;
    const int intersect = kmer_size > window_size;
	u64 intersect_mask = 0;
	int stride = kmer_size - window_size;
	if (intersect) intersect_mask = (U64_ONE << (stride << 1)) - 1;
    int kmer_idx = 0;

	if (!intersect) {
		for (size_t j = 0; j <= read_size - kmer_size; j += window_size) {
			u64 hash = 0;
			for (int k = 0; k < kmer_size; ++k) {
				size_t pos = j + k;
				u8 c = read[pos];
                hbn_assert(c >= 0 && c < 4);
				hash = (hash << 2) | c;
			}
            hash_array[kmer_idx++] = hash;
		}
	} else {
		u64 hash = 0;
		for (int j = 0; j < kmer_size; ++j) {
			size_t pos = j;
			u8 c = read[pos];
            hbn_assert(c >= 0 && c < 4);
			hash = (hash << 2) | c;
		}
		hash_array[kmer_idx++] = hash;
		for (u64 j = window_size; j <= read_size - kmer_size; j += window_size) {
			hash &= intersect_mask;
			for (int k = stride; k < kmer_size; ++k) {
				size_t pos = j + k;
                hbn_assert(pos < read_size, "p = %d, read_size = %d, j = %d, k = %d, stride = %d", pos, read_size, j, k, stride);
				u8 c = read[pos];
                hbn_assert(c >= 0 && c < 4);
				hash = (hash << 2) | c;
			}
			hash_array[kmer_idx++] = hash;
		}
	}
    return kmer_idx;
}

static int
collect_ddfkmer_subseq(const int context,
    const u8* read,
    const int read_from,
    const int read_to,
    const int kmer_size,
    const int window_size,
    vec_ddfk* ddfk_list)
{
    const int SL = kSeedingSeqSize, SR = kSeedingSeqStride;
    int n = read_to - read_from;
    int s = 0;
    u64 hash_list[kSeedingSeqSize];
    DDFKmer ddfk;
    int cnt = 0;
    ddfk.context = context;

    while (s < n) {
        int e = s + SL;
        e = hbn_min(e, n);
        int n_kmer = extract_hash_values(read + read_from + s, e - s, kmer_size, window_size, hash_list);
        for (int i = 0; i < n_kmer; ++i) {
            ddfk.offset = read_from + s + i * window_size;
            ddfk.hash = hash_list[i];
            kv_push(DDFKmer, *ddfk_list, ddfk);
            ++cnt;
        }
        s = e + SR;                
    }

    return cnt;
}

static int
collect_ddfkmer_one_query(const int context,
    const u8* read,
    const int read_size,
    const int kmer_size,
    const int kmer_window,
    vec_ddfk* ddfk_list)
{    
    return collect_ddfkmer_subseq(context, read, 0, read_size, kmer_size, kmer_window, ddfk_list);
}

static void
collect_ddfkmers(const text_t* query_vol,
    const int qidx_from,
    const int qidx_to,
    const int kmer_size,
    const int kmer_window,
    vec_ddfk* ddfk_list)
{
    kv_clear(*ddfk_list);
    kv_dinit(vec_u8, fwd_query);
    kv_dinit(vec_u8, rev_query);

    for (int i = qidx_from; i < qidx_to; ++i) {
        int ctx_id = i * 2;
        seqdb_extract_sequence(query_vol, i, FWD, &fwd_query);
        int query_length = kv_size(fwd_query);
        collect_ddfkmer_one_query(ctx_id, kv_data(fwd_query),
            query_length, kmer_size, kmer_window, ddfk_list);

        ++ctx_id;
        seqdb_extract_sequence(query_vol, i, REV, &rev_query);
        collect_ddfkmer_one_query(ctx_id, kv_data(rev_query),
            query_length, kmer_size, kmer_window, ddfk_list);
    }

    kv_destroy(fwd_query);
    kv_destroy(rev_query);
}

static int
proceed_to_next_ddfkm_idx(DDFKmer* ddfk_array, int ddfk_count, int ddfk_idx)
{
    if (ddfk_idx >= ddfk_count) return ddfk_idx;
    int hash = ddfk_array[ddfk_idx].hash;
    ++ddfk_idx;
    while (ddfk_idx < ddfk_count) {
        if (ddfk_array[ddfk_idx].hash != hash) break;
        ++ddfk_idx;
    }
    return ddfk_idx;
}

void*
collect_kmer_matches_thread(void* params)
{
    WordFinderThreadData* data = (WordFinderThreadData*)(params);
    kv_dinit(vec_ddfk, query_ddfk_list);
    collect_ddfkmers(data->query_vol, data->qidx_from, data->qidx_to, data->kmer_size, data->kmer_window, &query_ddfk_list);
    DDFKmer* query_ddfk_array = kv_data(query_ddfk_list);
    int query_ddfk_count = kv_size(query_ddfk_list);
    sort_ddfk_hash_lt(query_ddfk_count, query_ddfk_array);

    LookupTable* lktbl = data->lktbl;
    int qi = 0;
    int si = 0;
    kv_clear(data->l_ddfkm_list);
    while (qi < query_ddfk_count && si < lktbl->lki_count) {
        while (qi < query_ddfk_count && query_ddfk_array[qi].hash < lktbl->lki_array[si].hash) qi = proceed_to_next_ddfkm_idx(query_ddfk_array, query_ddfk_count, qi);
        if (qi >= query_ddfk_count) break;
        while (si < lktbl->lki_count && lktbl->lki_array[si].hash < query_ddfk_array[qi].hash) ++si;
        if (si >= lktbl->lki_count) break;
        if (query_ddfk_array[qi].hash == lktbl->lki_array[si].hash) {
            int next_qi = proceed_to_next_ddfkm_idx(query_ddfk_array, query_ddfk_count, qi);
            int qc = next_qi - qi;
            if (qc <= data->max_kmer_occ) {
                int sc = lktbl->lki_array[si].cnt;
                DDFKmer* qkm = query_ddfk_array + qi;
                u64* soff_array = lktbl->offset_list + lktbl->lki_array[si].offset_list_idx;
                DDFKmerMatch ddfkm;
                for (int i = 0; i < qc; ++i) {
                    ddfkm.context = qkm[i].context;
                    ddfkm.qoff = qkm[i].offset;
                    for (int j = 0; j < sc; ++j) {
                        ddfkm.soff = soff_array[j];
                        if (ddfkm.soff < data->global_soff_max_array[ddfkm.context]) {
                            kv_push(DDFKmerMatch, data->l_ddfkm_list, ddfkm);
                        }
                    }
                }
            }
            qi = next_qi;
            ++si;
        }
    }
    kv_destroy(query_ddfk_list);

    HBN_LOG("thread %d\t numberof kms: %d", data->thread_id, kv_size(data->l_ddfkm_list));

    DDFKmerMatch* km_array = kv_data(data->l_ddfkm_list);
    int km_count = kv_size(data->l_ddfkm_list);
    sort_ddfkm_context_lt(km_count, km_array);
    int i = 0;
    while (i < km_count) {
        int qid = km_array[i].context / 2;
        int max_ctx = qid * 2 + 1;
        int j = i + 1;
        while (j < km_count && km_array[j].context <= max_ctx) ++j;
        hbn_assert(qid < data->query_vol->dbinfo.num_seqs);
        data->g_kmi_array[qid].km_array = km_array + i;
        data->g_kmi_array[qid].km_count = j - i;
        i = j;
    }

    return NULL;
}

void
hbn_find_kmer_matches_mt(HbnWordFindData* data,
    const int qidx_from,
    const int qidx_to,
    const idx* soff_max_array)
{
    HbnWordFindData_FixThreadQidxRange(data, qidx_from, qidx_to);
    const int num_threads = data->num_threads;
    pthread_t jobid_array[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        data->wftd_array[i]->global_soff_max_array = soff_max_array;
        pthread_create(jobid_array + i, NULL, collect_kmer_matches_thread, data->wftd_array[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobid_array[i], NULL);
    }
}
