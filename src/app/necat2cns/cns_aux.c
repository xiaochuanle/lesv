#include "cns_aux.h"

#include "../../algo/hash_list_bucket_sort.h"
#include "../../corelib/partition_mt.h"
#include "../../ncbi_blast/c_ncbi_blast_aux.h"

CnsThreadData*
CnsThreadData_New(int thread_id,
    const HbnProgramOptions* opts,
    RawReadReader* raw_reads,
    RawReadCnsInfo* cns_info_array,
    int* cns_info_idx,
    pthread_mutex_t* cns_info_idx_lock,
    kstring_t* cns_fasta_array,
    pthread_mutex_t* out_lock)
{
    CnsThreadData* data = (CnsThreadData*)calloc(1, sizeof(CnsThreadData));
    data->thread_id = thread_id;
    data->opts = opts;
    data->raw_reads = raw_reads;
    data->cns_info_array = cns_info_array;
    data->cns_info_idx = cns_info_idx;
    data->cns_info_idx_lock = cns_info_idx_lock;
    ks_init(data->qaux);
    ks_init(data->saux);
    kv_init(data->fwd_read);
    kv_init(data->rev_read);
    kv_init(data->fwd_subject);
    kv_init(data->rev_subject);
    data->cns_fasta_array = cns_fasta_array;
    data->out_lock = out_lock;
    kv_init(data->cov_stats);
    data->cns_data = FCCnsDataNew();
    kv_init(data->chain_seed_list);
    data->hit_finder = InitHitFindDataNew(opts->memsc_kmer_size,
                            opts->memsc_kmer_window,
                            opts->memsc_score);
    data->tbck_data = HbnTracebackDataNew();
    data->ksw_data = Ksw2DataNew();
    ksw2_extd2_set_params(data->ksw_data);
    return data;
}

CnsThreadData*
CnsThreadData_Free(CnsThreadData* data)
{
    ks_destroy(data->qaux);
    ks_destroy(data->saux);
    kv_destroy(data->fwd_read);
    kv_destroy(data->rev_read);
    kv_destroy(data->fwd_subject);
    kv_destroy(data->rev_subject);
    kv_destroy(data->cov_stats);
    kv_destroy(data->chain_seed_list);
    data->cns_data = FCCnsDataFree(data->cns_data);
    data->hit_finder = InitHitFindDataFree(data->hit_finder);
    data->tbck_data = HbnTracebackDataFree(data->tbck_data);
    data->ksw_data = Ksw2DataFree(data->ksw_data);
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
load_and_sort_cns_hits(const char* can_dir, 
    const int pid, 
    const BOOL use_batch_mode,
    const int num_threads,
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
        if (cnt > MAX_CNS_SR) {
            ks_introsort_cns_hit_score_gt(cnt, hit_array + i);
            for (size_t k = i + MAX_CNS_SR; k < j; ++k) hit_array[k].sid = -1;
        }
        i = j;
    }

    i = 0;
    for (size_t k = 0; k < n; ++k) if (hit_array[k].sid >= 0) hit_array[i++] = hit_array[k];
    *hit_count = i;
    return hit_array;
}

BOOL
set_next_raw_read_batch_info(RawReadCnsInfo* cns_info_array,
    int* cns_info_count,
    HbnConsensusInitHit* cns_hit_array,
    size_t cns_hit_count,
    size_t* cns_hit_idx,
    RawReadReader* raw_reads,
    const int batch_size)
{
    size_t i = *cns_hit_idx;
    if (i >= cns_hit_count) return FALSE;
    int p = 0;
    memset(cns_info_array, 0, sizeof(RawReadCnsInfo) * batch_size);
    int hit_cnt = 0;
    while (i < cns_hit_count) {
        size_t j = i + 1;
        while (j < cns_hit_count && cns_hit_array[i].sid == cns_hit_array[j].sid) ++j;
        cns_info_array[p].oid = cns_hit_array[i].sid;
        cns_info_array[p].can_from = i;
        cns_info_array[p].can_to = j;
        hit_cnt += j - i;
        i = j;
        ++p;
        if (p == batch_size) break;
    }

    *cns_info_count = p;
    *cns_hit_idx = i;
    RawReadReader_LoadRawReadFromCnsHitArray(cns_hit_array + cns_info_array[0].can_from,
        hit_cnt,
        raw_reads);
    return TRUE;
}