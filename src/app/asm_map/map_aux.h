#ifndef __MAP_AUX_H
#define __MAP_AUX_H

#include "hbn_options.h"
#include "../../algo/init_hit_finder.h"
#include "../../algo/hbn_traceback.h"
#include "../../corelib/raw_reads_reader.h"
#include "../../ncbi_blast/setup/blast_hits.h"

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    BlastHSPList* hsplist_array;
    int hsplist_count;
    BlastHSP** hsp_ptr_array;
    BlastHSP* hsp_array;
    int max_sr;
    int max_hsp_perc_sr;
    kstring_t align_string;
} AsmMapResults;

void
AsmMapResultsClear(AsmMapResults* results);

AsmMapResults*
AsmMapResultsNew(int max_sr, int max_hsp_perc_sr);

AsmMapResults*
AsmMapResultsFree(AsmMapResults* results);

typedef struct {
    int oid;
    size_t sr_idx_from;
    size_t sr_idx_to;
} SubjectSupportReadInfo;

typedef struct {
    int thread_id;
    const HbnProgramOptions* opts;
    RawReadReader* raw_reads;
    SubjectSupportReadInfo* sr_info_array;
    int sr_info_count;
    int* sr_info_idx;
    pthread_mutex_t* sr_info_idx_lock;
    HbnConsensusInitHit* hit_array;
    size_t hit_count;
    kstring_t qaux;
    kstring_t saux;
    vec_u8 fwd_read;
    vec_u8 rev_read;
    vec_u8 fwd_subject;
    vec_u8 rev_subject;
    kstring_t out_buf;
    FILE* out;
    pthread_mutex_t* out_lock;
    vec_chain_seed chain_seed_list;
    InitHitFindData* hit_finder;
    HbnTracebackData* tbck_data;
    Ksw2Data* ksw_data;
    AsmMapResults* results;
} AsmMapThreadData;

AsmMapThreadData*
AsmMapThreadDataNew(int thread_id,
    const HbnProgramOptions* opts,
    RawReadReader* raw_reads,
    SubjectSupportReadInfo* sr_info_array,
    int* sr_info_idx,
    pthread_mutex_t* sr_info_idx_lock,
    FILE* out,
    pthread_mutex_t* out_lock);

AsmMapThreadData*
AsmMapThreadDataFree(AsmMapThreadData* data);

BOOL
set_next_subject_batch_info(SubjectSupportReadInfo* sr_info_array,
    int* sr_info_count,
    HbnConsensusInitHit* sr_hit_array,
    size_t sr_hit_count,
    size_t* sr_hit_idx,
    RawReadReader* raw_reads,
    const int batch_size);

HbnConsensusInitHit*
load_and_sort_sr_hits(const char* can_dir, 
    const int pid, 
    const BOOL use_batch_mode,
    const int num_threads,
    const int max_sr_per_subject,
    size_t* hit_count);

#ifdef __cplusplus
}
#endif

#endif // __MAP_AUX_H