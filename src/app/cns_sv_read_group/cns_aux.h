#ifndef __CNS_AUX_H
#define __CNS_AUX_H

#include "../../algo/fccns/fccns.h"
#include "../../algo/hbn_traceback.h"
#include "../../algo/init_hit_finder.h"
#include "../../corelib/seqdb.h"
#include "cns_options.h"
#include "sv_read_group.h"

#ifdef __cplusplus
extern "C" {
#endif

#define INDEL_COV_FACTOR (0.4)

typedef struct {
    HbnProgramOptions* opts;
    CSeqDB* db;
    const CSeqInfo* raw_read_info_array;
    const char* raw_read_name_array;
    InitHitFindData* hit_finder;
    HbnTracebackData* tbck_data;
    Ksw2Data* ksw_data;
    FCCnsData* cns_data;
    FILE* fasta_out;
    FILE* sam_out;
    pthread_mutex_t* out_lock;
    int* next_svr_id;
    pthread_mutex_t* svr_id_lock;
    kstring_t qaln;
    kstring_t saln;
    vec_chain_seed chain_seed_list;
    SvReadGroup* raw_group;
    SvReadGroup* cns_group;
    int cns_iter;
} CnsThreadData;

CnsThreadData*
CnsThreadDataNew(HbnProgramOptions* opts, 
    CSeqDB* db,
    const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    FILE* fasta_out,
    FILE* sam_out,
    pthread_mutex_t* out_lock,
    int* svr_id,
    pthread_mutex_t* svr_id_lock);

CnsThreadData*
CnsThreadDataFree(CnsThreadData* data);

#ifdef __cplusplus
}
#endif

#endif // __CNS_AUX_H