#include "correct_one_part.h"

#include "../../corelib/partition_mt.h"
#include "../../ncbi_blast/c_ncbi_blast_aux.h"
#include "correct_one_read.h"

typedef struct {
    int pid;
    const HbnProgramOptions* opts;
    RawReadReader* raw_reads;
    HbnConsensusInitHit* hit_array;
    size_t hit_count;
    size_t hit_idx;
    RawReadCnsInfo* cns_info_array;
    int cns_info_count;
    int cns_info_idx;
    pthread_mutex_t cns_info_idx_lock;
    kstring_t cns_fasta_array;
    pthread_mutex_t cns_fasta_array_lock;
    FILE* out;
    CnsThreadData** data_array;
    FILE* cns_seq_id_out;
} CnsOnePartData;

static CnsOnePartData*
CnsOnePartDataNew(const HbnProgramOptions* opts,
    RawReadReader* raw_reads,
    int pid)
{
    CnsOnePartData* data = (CnsOnePartData*)calloc(1, sizeof(CnsOnePartData));
    data->pid = pid;
    data->opts = opts;
    data->raw_reads = raw_reads;
    data->cns_info_array = (RawReadCnsInfo*)calloc(opts->batch_size, sizeof(RawReadCnsInfo));
    data->hit_array = load_and_sort_cns_hits(opts->can_dir, 
                        pid, 
                        opts->use_batch_mode, 
                        opts->num_threads, 
                        &data->hit_count);
    for (int i = 0; i < data->hit_count; ++i) {
        hbn_assert(data->hit_array[i].qid < raw_reads->dbinfo.num_seqs);
        hbn_assert(data->hit_array[i].sid < raw_reads->dbinfo.num_seqs);
    }
    HBN_LOG("load %zu hits", data->hit_count);
    pthread_mutex_init(&data->cns_info_idx_lock, NULL);
    ks_init(data->cns_fasta_array);
    pthread_mutex_init(&data->cns_fasta_array_lock, NULL);
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(opts->can_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".cns.fasta");
    HBN_LOG("corrected reads path: %s", path);
    hbn_fopen(data->out, path, "w");
    make_partition_name(opts->can_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".cns_seq_id");
    hbn_fopen(data->cns_seq_id_out, path, "w");
    data->data_array = (CnsThreadData**)calloc(opts->num_threads, sizeof(CnsThreadData*));
    for (int i = 0; i < opts->num_threads; ++i) {
        data->data_array[i] = CnsThreadData_New(i,
                                opts,
                                raw_reads,
                                data->cns_info_array,
                                &data->cns_info_idx,
                                &data->cns_info_idx_lock,
                                &data->cns_fasta_array,
                                &data->cns_fasta_array_lock);
        data->data_array[i]->hit_array = data->hit_array;
        data->data_array[i]->hit_count = data->hit_count;
    }
    return data;
}

static CnsOnePartData*
CnsOnePartDataFree(CnsOnePartData* data)
{
    if (!data) return NULL;
    if (data->hit_array) sfree(data->hit_array);
    data->hit_count = 0;
    if (data->cns_info_array) sfree(data->cns_info_array);
    ks_destroy(data->cns_fasta_array);
    if (data->out) hbn_fclose(data->out);
    data->out = NULL;
    if (data->cns_seq_id_out) hbn_fclose(data->cns_seq_id_out);
    data->cns_seq_id_out = NULL;
    for (int i = 0; i < data->opts->num_threads; ++i) {
        data->data_array[i] = CnsThreadData_Free(data->data_array[i]);
    }
    sfree(data->data_array);
    sfree(data);
    return NULL;
}

static BOOL
CnsOnePartData_SetupNextCanBatch(CnsOnePartData* data)
{
    if (!set_next_raw_read_batch_info(data->cns_info_array,
            &data->cns_info_count,
            data->hit_array,
            data->hit_count,
            &data->hit_idx,
            data->raw_reads,
            data->opts->batch_size)) {
        return FALSE;
    }    
    data->cns_info_idx = 0;
    ks_clear(data->cns_fasta_array);
    for (int i = 0; i < data->opts->num_threads; ++i) {
        data->data_array[i]->cns_info_count = data->cns_info_count;
    }
    return TRUE;
}

static void
CnsOnePartData_DumpCnsFasta(CnsOnePartData* data)
{
    for (int i = 0; i < data->cns_info_count; ++i) {
        RawReadCnsInfo* cns_info = data->cns_info_array + i;
        if (cns_info->cns_fasta_size == 0) {
            continue;
        } else {
            const char* fasta = ks_s(data->cns_fasta_array) + cns_info->cns_fasta_offset;
            hbn_fwrite(fasta, 1, cns_info->cns_fasta_size, data->out);
            hbn_fwrite(&cns_info->oid, sizeof(int), 1, data->cns_seq_id_out);
        }
    }
}

static BOOL
partition_is_corrected(const char* can_dir, const int pid)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(can_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".corrected");
    return access(path, F_OK) == 0;
}

static void
partition_make_corrected(const char* can_dir, const int pid)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(can_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".corrected");
    hbn_dfopen(out, path, "w");
    hbn_fclose(out);
}

static void*
cns_thread_worker(void* params)
{
    CnsThreadData* data = (CnsThreadData*)(params);
    while (1) {
        int cns_info_idx = -1;
        pthread_mutex_lock(data->cns_info_idx_lock);
        cns_info_idx = *data->cns_info_idx;
        ++(*data->cns_info_idx);
        pthread_mutex_unlock(data->cns_info_idx_lock);
        if (cns_info_idx >= data->cns_info_count) break;
        correct_one_read(data, cns_info_idx);
    }
    return NULL;
}

void
correct_one_part(const HbnProgramOptions* opts, RawReadReader* raw_reads, const int pid)
{
    if (partition_is_corrected(opts->can_dir, pid)) {
        HBN_LOG("Partition %d has been corrected. Skip it.");
        return;
    }
    char job_name[512];
    char sid1_str[64];
    char sid2_str[64];
    pthread_t jobs[opts->num_threads];
    CnsOnePartData* data = CnsOnePartDataNew(opts, raw_reads, pid);
    while (CnsOnePartData_SetupNextCanBatch(data)) {
        int sid1 = data->cns_info_array[0].oid;
        int sid2 = data->cns_info_array[data->cns_info_count-1].oid + 1;
	//if (sid1 < 36000) continue;
        u64_to_fixed_width_string_r(sid1, sid1_str, HBN_DIGIT_WIDTH);
        u64_to_fixed_width_string_r(sid2, sid2_str, HBN_DIGIT_WIDTH);
        sprintf(job_name, "Correct subject %s --- %s", sid1_str, sid2_str);
        hbn_timing_begin(job_name);
        for (int i = 0; i < opts->num_threads; ++i) {
            pthread_create(jobs + i, NULL, cns_thread_worker, data->data_array[i]);
        }
        for (int i = 0; i < opts->num_threads; ++i) {
            pthread_join(jobs[i], NULL);
        }
        CnsOnePartData_DumpCnsFasta(data);
        hbn_timing_end(job_name);
    }
    data = CnsOnePartDataFree(data);
    partition_make_corrected(opts->can_dir, pid);
}
