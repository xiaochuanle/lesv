#include "map_aux.h"

#include "../../algo/hash_list_bucket_sort.h"
#include "../../corelib/partition_mt.h"
#include "../../ncbi_blast/c_ncbi_blast_aux.h"

void
AsmMapResultsClear(AsmMapResults* results)
{
    for (int i = 0; i < results->hsplist_count; ++i) {
        results->hsplist_array[i].hspcnt = 0;
        results->hsplist_array[i].oid = -1;
        results->hsplist_array[i].hbn_best_raw_score = 0;
    }
    ks_clear(results->align_string);
}

AsmMapResults*
AsmMapResultsNew(int max_sr, int max_hsp_perc_sr)
{
    AsmMapResults* results = (AsmMapResults*)calloc(1, sizeof(AsmMapResults));
    results->hsplist_array = (BlastHSPList*)calloc(max_sr, sizeof(BlastHSPList));
    results->hsplist_count = 0;
    results->hsp_ptr_array = (BlastHSP**)calloc(max_sr * max_hsp_perc_sr, sizeof(BlastHSP*));
    results->hsp_array = (BlastHSP*)calloc(max_sr * max_hsp_perc_sr, sizeof(BlastHSP));
    results->max_sr = max_sr;
    results->max_hsp_perc_sr = max_hsp_perc_sr;
    ks_init(results->align_string);

    for (int i = 0; i < max_sr * max_hsp_perc_sr; ++i) results->hsp_ptr_array[i] = &results->hsp_array[i];

    BlastHSP** p = results->hsp_ptr_array;
    for (int i = 0; i < max_sr; ++i) {
        results->hsplist_array[i].oid = -1;
        results->hsplist_array[i].hspcnt = 0;
        results->hsplist_array[i].hsp_max = max_hsp_perc_sr;
        results->hsplist_array[i].hbn_best_raw_score = 0;
        results->hsplist_array[i].hsp_array = p;
        p += max_hsp_perc_sr;
    }

    return results;
}

AsmMapResults*
AsmMapResultsFree(AsmMapResults* results)
{
    free(results->hsplist_array);
    free(results->hsp_ptr_array);
    free(results->hsp_array);
    ks_destroy(results->align_string);
    free(results);
    return NULL;
}

AsmMapThreadData*
AsmMapThreadDataNew(int thread_id,
    const HbnProgramOptions* opts,
    RawReadReader* raw_reads,
    SubjectSupportReadInfo* sr_info_array,
    int* sr_info_idx,
    pthread_mutex_t* sr_info_idx_lock,
    FILE* out,
    pthread_mutex_t* out_lock)
{
    AsmMapThreadData* data = (AsmMapThreadData*)calloc(1, sizeof(AsmMapThreadData));
    data->thread_id = thread_id;
    data->opts = opts;
    data->raw_reads = raw_reads;
    data->sr_info_array = sr_info_array;
    data->sr_info_idx = sr_info_idx;
    data->sr_info_idx_lock = sr_info_idx_lock;
    ks_init(data->qaux);
    ks_init(data->saux);
    kv_init(data->fwd_read);
    kv_init(data->rev_read);
    kv_init(data->fwd_subject);
    kv_init(data->rev_subject);
    ks_init(data->out_buf);
    data->out = out;
    data->out_lock = out_lock;
    kv_init(data->chain_seed_list);
    data->hit_finder = InitHitFindDataNew(opts->memsc_kmer_size,
                            opts->memsc_kmer_window,
                            opts->memsc_score);
    data->tbck_data = HbnTracebackDataNew();
    data->ksw_data = Ksw2DataNew();
    ksw2_extd2_set_params(data->ksw_data);
    data->results = AsmMapResultsNew(opts->hitlist_size, opts->max_hsps_per_subject);
    return data;
}

AsmMapThreadData*
AsmMapThreadDataFree(AsmMapThreadData* data)
{
    ks_destroy(data->qaux);
    ks_destroy(data->saux);
    kv_destroy(data->fwd_read);
    kv_destroy(data->rev_read);
    kv_destroy(data->fwd_subject);
    kv_destroy(data->rev_subject);
    ks_destroy(data->out_buf);
    kv_destroy(data->chain_seed_list);
    data->hit_finder = InitHitFindDataFree(data->hit_finder);
    data->tbck_data = HbnTracebackDataFree(data->tbck_data);
    data->ksw_data = Ksw2DataFree(data->ksw_data);
    data->results = AsmMapResultsFree(data->results);
    sfree(data);
    return data;
}

//// radix sort related functions

u64 radix_sort_cns_hit_sid(void* a, const u64 i)
{
    HbnConsensusInitHit* p = (HbnConsensusInitHit*)(a);
    return p[i].sid;
}

static void
set_cns_hit_array_item_value(void* src, const u64 src_idx, void* dst, const u64 dst_idx)
{
    HbnConsensusInitHit* src_array = (HbnConsensusInitHit*)(src);
    HbnConsensusInitHit* dst_array = (HbnConsensusInitHit*)(dst);
    dst_array[dst_idx] = src_array[src_idx];
}

static void
sort_cns_hit_array_sid_lt_mt(size_t n, void* a, int num_threads)
{
    radix_sort(a, 
		   sizeof(HbnConsensusInitHit),
		   n, 
		   num_threads,
		   NULL, 
		   radix_sort_cns_hit_sid,
		   set_cns_hit_array_item_value);
}

/// end radix sort related functions

HbnConsensusInitHit*
load_and_sort_sr_hits(const char* can_dir, 
    const int pid, 
    const BOOL use_batch_mode,
    const int num_threads,
    const int max_sr_per_subject,
    size_t* hit_count)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(can_dir, DEFAULT_PART_PREFIX, pid, path);
    HBN_LOG("can_path: %s", path);
    HbnConsensusInitHit* hit_array = (HbnConsensusInitHit*)load_part_records(path, sizeof(HbnConsensusInitHit), hit_count);
    size_t n = *hit_count;
    if (n == 0) return NULL;
    sort_cns_hit_array_sid_lt_mt(n, hit_array, num_threads);
    if (!use_batch_mode) return hit_array;

    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && hit_array[i].sid == hit_array[j].sid) ++j;
        int cnt = j - i;
        if (cnt > max_sr_per_subject) {
            ks_introsort_cns_hit_score_gt(cnt, hit_array + i);
            for (size_t k = i + max_sr_per_subject; k < j; ++k) hit_array[k].sid = -1;
        }
        i = j;
    }

    i = 0;
    for (size_t k = 0; k < n; ++k) if (hit_array[k].sid >= 0) hit_array[i++] = hit_array[k];
    *hit_count = i;
    return hit_array;
}

BOOL
set_next_subject_batch_info(SubjectSupportReadInfo* sr_info_array,
    int* sr_info_count,
    HbnConsensusInitHit* sr_hit_array,
    size_t sr_hit_count,
    size_t* sr_hit_idx,
    RawReadReader* raw_reads,
    const int batch_size)
{
    size_t i = *sr_hit_idx;
    if (i >= sr_hit_count) return FALSE;
    int p = 0;
    memset(sr_info_array, 0, sizeof(SubjectSupportReadInfo) * batch_size);
    int hit_cnt = 0;
    while (i < sr_hit_count) {
        size_t j = i + 1;
        while (j < sr_hit_count && sr_hit_array[i].sid == sr_hit_array[j].sid) ++j;
        sr_info_array[p].oid = sr_hit_array[i].sid;
        sr_info_array[p].sr_idx_from = i;
        sr_info_array[p].sr_idx_to = j;
        hit_cnt += j - i;
        i = j;
        ++p;
        if (p == batch_size) break;
    }

    *sr_info_count = p;
    *sr_hit_idx = i;
    RawReadReader_LoadRawReadFromCnsHitArray(sr_hit_array + sr_info_array[0].sr_idx_from,
        hit_cnt,
        raw_reads);
    return TRUE;
}