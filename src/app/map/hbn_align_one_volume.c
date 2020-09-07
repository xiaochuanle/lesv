#include "hbn_align_one_volume.h"

#include "chain_and_extend_kmer_matches.h"
#include "hbn_find_subseq_hit.h"
#include "hbn_extend_subseq_hit.h"
#include "hbn_extend_subseq_hit_diff.h"
#include "mecat_results.h"
#include "../../corelib/m4_record.h"
#include "../../corelib/collect_km_mode.h"
#include "../../ncbi_blast/setup/hsp2string.h"

#include <pthread.h>

static hbn_task_struct* g_task_struct = NULL;
static int g_thread_index = 0;
static pthread_mutex_t g_thread_index_lock;
static int g_query_index = 0;
static pthread_mutex_t g_query_index_lock;
static int g_thread_id = 0;
static pthread_mutex_t g_thread_id_lock;

static size_t gBatchQueryRes = 0;
static size_t g_left_query_res = 0;
static int g_batch_qidx_from = 0;
static int g_batch_qidx_to = 0;

static int g_max_query_global_id = -1;
static int g_max_subject_global_id = -1;

static void
init_global_values(hbn_task_struct* task_struct)
{
    g_task_struct = task_struct;
    g_thread_index = 0;
    pthread_mutex_init(&g_thread_index_lock, NULL);
    g_query_index = 0;
    pthread_mutex_init(&g_query_index_lock, NULL);
    g_thread_id = 0;
    pthread_mutex_init(&g_thread_id_lock, NULL);
    g_batch_qidx_from = 0;
    g_batch_qidx_to = 0;
    for (int i = 0; i < g_task_struct->query_vol->dbinfo.num_seqs; ++i) {
        g_left_query_res += seqdb_seq_size(g_task_struct->query_vol, i);
    }
    gBatchQueryRes = task_struct->opts->query_batch_size;

    g_max_query_global_id = seqdb_load_num_reads(task_struct->opts->db_dir, INIT_QUERY_DB_TITLE);
    if (task_struct->opts->align_task == eHbnTask_pm) {
        g_max_subject_global_id = g_max_query_global_id;
    } else {
        g_max_subject_global_id = seqdb_load_num_reads(task_struct->opts->db_dir, INIT_SUBJECT_DB_TITLE);
    }
}

static BOOL
get_next_batch_query_qidx_range()
{
    const text_t* query_vol = g_task_struct->query_vol;
    int nq = seqdb_num_seqs(query_vol);
    int from = g_batch_qidx_to;
    if (from >= nq) return FALSE;
    size_t curr = seqdb_seq_size(query_vol, from);
    int to = from + 1;
    while (curr < gBatchQueryRes && to < nq) {
        curr += seqdb_seq_size(query_vol, to);
        ++to;
    }

    size_t kMinBatchRes = gBatchQueryRes * 0.2;
    hbn_assert(g_left_query_res >= curr);
    g_left_query_res -= curr;
    if (g_left_query_res <= kMinBatchRes) {
        to = nq;
        g_left_query_res = 0;
    }

    char res_buf[64];
    u64_to_string_datasize(curr, res_buf);
    HBN_LOG("mapping queries %d --- %d (%d, %s)", from, to, to - from, res_buf);
    g_batch_qidx_from = from;
    g_batch_qidx_to = to;
    return TRUE;
}

static idx*
fill_max_soff_array(const text_t* query_vol,
    const text_t* subject_vol,
    const BOOL query_and_subject_are_the_same)
{
    int num_contexts = 2 * query_vol->dbinfo.num_seqs;
    idx* soff_max_array = (idx*)malloc(sizeof(idx) * num_contexts);
    for (int i = 0; i < num_contexts; ++i) soff_max_array[i] = IDX_MAX;
    if (!query_and_subject_are_the_same) return soff_max_array;

    const int reference_start_id = subject_vol->dbinfo.seq_start_id;
    const int read_start_id = query_vol->dbinfo.seq_start_id;
    int max_refid = reference_start_id + seqdb_num_seqs(subject_vol);
    for (int i = 0; i < query_vol->dbinfo.num_seqs; ++i) {
        int read_id = i;
        int g_read_id = read_id + read_start_id;
        if (g_read_id >= reference_start_id && g_read_id < max_refid) {
            idx soff_max = seqdb_seq_offset(query_vol, read_id);
            soff_max_array[2*i] = soff_max;
            soff_max_array[2*i+1] = soff_max;
        }
    }
    return soff_max_array;
}

void
hbn_align_one_volume(hbn_task_struct* task_struct)
{
    init_global_values(task_struct);
    idx* soff_max_array = fill_max_soff_array(task_struct->query_vol,
                            task_struct->subject_vol,
                            task_struct->query_and_subject_are_the_same);
    
    char job_name[256];
    char qfstr[64], qtstr[64];
    while (get_next_batch_query_qidx_range()) {
        u64_to_fixed_width_string_r(g_batch_qidx_from, qfstr, HBN_DIGIT_WIDTH);
        u64_to_fixed_width_string_r(g_batch_qidx_to, qtstr, HBN_DIGIT_WIDTH);
        sprintf(job_name, "mapping query %s --- %s", qfstr, qtstr);
        hbn_timing_begin(job_name);
        HbnWordFindData* word_data = HbnWordFindDataNew(task_struct->lktbl,
                                        task_struct->subject_vol,
                                        task_struct->query_vol,
                                        task_struct->opts->num_threads,
                                        task_struct->opts->kmer_size,
                                        1,
                                        task_struct->opts->max_kmer_occ);
        hbn_find_kmer_matches_mt(word_data, g_batch_qidx_from, g_batch_qidx_to, soff_max_array);
        chain_and_extend_km_mt(word_data, task_struct->opts, 
            task_struct->query_vol,
            task_struct->subject_vol,
            g_batch_qidx_from,
            g_batch_qidx_to,
            task_struct->out,
            task_struct->qi_vs_sj_out);
        word_data = HbnWordFindDataFree(word_data);   
        hbn_timing_end(job_name);
    }                       
    free(soff_max_array);
}
