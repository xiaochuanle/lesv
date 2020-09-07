#ifndef __CNS_AUX_H
#define __CNS_AUX_H

#include "cns_options.h"
#include "../../algo/fccns/fccns.h"
#include "../../algo/hbn_traceback.h"
#include "../../algo/init_hit_finder.h"
#include "../../corelib/raw_reads_reader.h"

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_CNS_SR  80
#define MAX_CNS_COV 12

#define INDEL_COV_FACTOR (0.1)

typedef struct {
    int oid;
    int raw_read_size;
    int cns_from;
    int cns_to;
    int cns_read_size;
    size_t can_from;
    size_t can_to;
    size_t cns_fasta_offset;
    size_t cns_fasta_size;
} RawReadCnsInfo;

typedef struct {
    int thread_id;
    const HbnProgramOptions* opts;
    RawReadReader* raw_reads;
    RawReadCnsInfo* cns_info_array;
    int cns_info_count;
    int* cns_info_idx;
    pthread_mutex_t* cns_info_idx_lock;
    HbnConsensusInitHit* hit_array;
    size_t hit_count;
    kstring_t qaux;
    kstring_t saux;
    vec_u8 fwd_read;
    vec_u8 rev_read;
    vec_u8 fwd_subject;
    vec_u8 rev_subject;
    kstring_t* cns_fasta_array;
    pthread_mutex_t* out_lock;
    vec_u8 cov_stats;
    FCCnsData* cns_data;
    vec_chain_seed chain_seed_list;
    InitHitFindData* hit_finder;
    HbnTracebackData* tbck_data;
    Ksw2Data* ksw_data;
} CnsThreadData;

CnsThreadData*
CnsThreadData_New(int thread_id,
    const HbnProgramOptions* opts,
    RawReadReader* raw_reads,
    RawReadCnsInfo* cns_info_array,
    int* cns_info_idx,
    pthread_mutex_t* cns_info_idx_lock,
    kstring_t* cns_fasta_array,
    pthread_mutex_t* out_lock);

CnsThreadData*
CnsThreadData_Free(CnsThreadData* data);

HbnConsensusInitHit*
load_and_sort_cns_hits(const char* can_dir, 
    const int pid, 
    const BOOL use_batch_mode,
    const int num_threads,
    size_t* hit_count);

BOOL
set_next_raw_read_batch_info(RawReadCnsInfo* cns_info_array,
    int* cns_info_count,
    HbnConsensusInitHit* cns_hit_array,
    size_t cns_hit_count,
    size_t* cns_hit_idx,
    RawReadReader* raw_reads,
    const int batch_size);

#ifdef __cplusplus
}
#endif

#endif // __CNS_AUX_H
