#include "map_one_part.h"

#include "../../corelib/partition_mt.h"
#include "map_one_read.h"

typedef struct {
    int pid;
    const HbnProgramOptions* opts;
    RawReadReader* raw_reads;
    HbnConsensusInitHit* hit_array;
    size_t hit_count;
    size_t hit_idx;
    int batch_size;
    SubjectSupportReadInfo* sr_info_array;
    int sr_info_count;
    int sr_info_idx;
    pthread_mutex_t sr_info_idx_lock;
    FILE* out;
    pthread_mutex_t out_lock;
    AsmMapThreadData** data_array;
} AsmMapOnePartData;

static AsmMapOnePartData*
AsmMapOnePartDataNew(const char* psr_dir,
    const HbnProgramOptions* opts,
    RawReadReader* raw_reads,
    const int pid,
    const int batch_size,
    const BOOL use_batch_mode)
{
    AsmMapOnePartData* data = (AsmMapOnePartData*)calloc(1, sizeof(AsmMapOnePartData));
    data->pid = pid;
    data->opts = opts;
    data->raw_reads = raw_reads;
    data->batch_size = batch_size;
    data->sr_info_array = (SubjectSupportReadInfo*)calloc(batch_size, sizeof(SubjectSupportReadInfo));
    data->hit_array = load_and_sort_sr_hits(psr_dir, 
                        pid, 
                        use_batch_mode, 
                        opts->num_threads, 
                        opts->hitlist_size,
                        &data->hit_count);
    for (int i = 0; i < data->hit_count; ++i) {
        hbn_assert(data->hit_array[i].qid < raw_reads->dbinfo.num_seqs);
        hbn_assert(data->hit_array[i].sid < raw_reads->dbinfo.num_seqs);
    }
    HBN_LOG("load %zu hits", data->hit_count);
    pthread_mutex_init(&data->sr_info_idx_lock, NULL);
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(psr_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".ovlp");
    hbn_fopen(data->out, path, "w");
    pthread_mutex_init(&data->out_lock, NULL);
    data->data_array = (AsmMapThreadData**)calloc(opts->num_threads, sizeof(AsmMapThreadData*));
    for (int i = 0; i < opts->num_threads; ++i) {
        data->data_array[i] = AsmMapThreadDataNew(i,
                                opts,
                                raw_reads,
                                data->sr_info_array,
                                &data->sr_info_idx,
                                &data->sr_info_idx_lock,
                                data->out,
                                &data->out_lock);
        data->data_array[i]->hit_array = data->hit_array;
        data->data_array[i]->hit_count = data->hit_count;
    }
    return data;
}

static AsmMapOnePartData*
AsmMapOnePartDataFree(AsmMapOnePartData* data)
{
    if (!data) return NULL;
    if (data->hit_array) sfree(data->hit_array);
    data->hit_count = 0;
    if (data->sr_info_array) sfree(data->sr_info_array);
    if (data->out) hbn_fclose(data->out);
    data->out = NULL;
    for (int i = 0; i < data->opts->num_threads; ++i) {
        data->data_array[i] = AsmMapThreadDataFree(data->data_array[i]);
    }
    sfree(data->data_array);
    sfree(data);
    return NULL;
}

static BOOL
AsmMapOnePartData_SetupNextSrBatch(AsmMapOnePartData* data)
{
    if (!set_next_subject_batch_info(data->sr_info_array,
            &data->sr_info_count,
            data->hit_array,
            data->hit_count,
            &data->hit_idx,
            data->raw_reads,
            data->batch_size)) {
        return FALSE;
    }    
    data->sr_info_idx = 0;
    for (int i = 0; i < data->opts->num_threads; ++i) {
        data->data_array[i]->sr_info_count = data->sr_info_count;
    }
    return TRUE;
}

static BOOL
partition_is_mapped(const char* can_dir, const int pid)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(can_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".mapped");
    return access(path, F_OK) == 0;
}

static void
partition_make_mapped(const char* can_dir, const int pid)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(can_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".mapped");
    hbn_dfopen(out, path, "w");
    hbn_fclose(out);
}

static void*
asm_map_thread_worker(void* params)
{
    AsmMapThreadData* data = (AsmMapThreadData*)(params);
    while (1) {
        int sr_info_idx = -1;
        pthread_mutex_lock(data->sr_info_idx_lock);
        sr_info_idx = *data->sr_info_idx;
        ++(*data->sr_info_idx);
        pthread_mutex_unlock(data->sr_info_idx_lock);
        if (sr_info_idx >= data->sr_info_count) break;
        asm_map_one_read(data, sr_info_idx);
    }
    return NULL;
}

void
asm_map_one_part(const char* psr_dir,
    const HbnProgramOptions* opts, 
    RawReadReader* raw_reads, 
    const int pid,
    const int batch_size,
    const BOOL use_batch_mode)
{
    if (partition_is_mapped(psr_dir, pid)) {
        HBN_LOG("Partition %d has been mapped. Skip it.");
        return;
    }
    char job_name[512];
    char sid1_str[64];
    char sid2_str[64];
    pthread_t jobs[opts->num_threads];
    AsmMapOnePartData* data = AsmMapOnePartDataNew(psr_dir, opts, raw_reads, pid, batch_size, use_batch_mode);
    while (AsmMapOnePartData_SetupNextSrBatch(data)) {
        int sid1 = data->sr_info_array[0].oid;
        int sid2 = data->sr_info_array[data->sr_info_count-1].oid + 1;
        u64_to_fixed_width_string_r(sid1, sid1_str, HBN_DIGIT_WIDTH);
        u64_to_fixed_width_string_r(sid2, sid2_str, HBN_DIGIT_WIDTH);
        sprintf(job_name, "Map subject %s --- %s", sid1_str, sid2_str);
        hbn_timing_begin(job_name);
        for (int i = 0; i < opts->num_threads; ++i) {
            pthread_create(jobs + i, NULL, asm_map_thread_worker, data->data_array[i]);
        }
        for (int i = 0; i < opts->num_threads; ++i) {
            pthread_join(jobs[i], NULL);
        }
        hbn_timing_end(job_name);
    }
    data = AsmMapOnePartDataFree(data);
    partition_make_mapped(psr_dir, pid);
}  