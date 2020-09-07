#ifndef __CNS_EXTEND_CHAIN_SEED_LIST_H
#define __CNS_EXTEND_CHAIN_SEED_LIST_H

#include "hbn_traceback.h"
#include "hbn_traceback_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

BOOL
check_ovlp_position(const int qb, const int qe, const int qs,
    const int sb, const int se, const int ss);

BOOL
check_ovlp_cov(const int qb, const int qe, const int qs,
						const int sb, const int se, const int ss,
                        int res_coverage,
						double perc_coverage);

BOOL
check_chain_seed_array_cov(ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const int read_length,
    const int subject_length,
    const int res_coverage,
    const  double perc_coverage);

BOOL
subject_subseq_cov_is_full(u8* cov_stats, int soff, int send, const int kMaxCov);

BOOL
cns_extend_chain_seed_array(HbnTracebackData* traceback_data,
    Ksw2Data* ksw_data,
    ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const int min_cov_res,
    const double min_perc_coverage,
    const double min_perc_identity,
    const int read_id,
    const u8* read,
    const int read_length,
    const int subject_id,
    const u8* subject,
    const int subject_length,
    kstring_t* qaln,
    kstring_t* saln,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* ident_perc_);

BOOL
cns_extend_chain_seed_array_with_refine(HbnTracebackData* traceback_data,
    Ksw2Data* ksw_data,
    ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const int min_cov_res,
    const double min_perc_coverage,
    const double min_perc_identity,
    const int read_id,
    const u8* read,
    const int read_length,
    const int subject_id,
    const u8* subject,
    const int subject_length,
    kstring_t* qaln,
    kstring_t* saln,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* ident_perc_);

BOOL 
cns_extend_chain_seed_array_with_ksw(HbnTracebackData* traceback_data,
    Ksw2Data* ksw2_data,
    ChainSeed* chain_seed_array,
    const int chain_seed_count,
    const int min_cov_res,
    const double min_perc_coverage,
    const double min_perc_identity,
    const int read_id,
    const u8* read,
    const int read_length,
    const int subject_id,
    const u8* subject,
    const int subject_length,
    kstring_t* qaln,
    kstring_t* saln,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* ident_perc_);

#ifdef __cplusplus
}
#endif

#endif // __CNS_EXTEND_CHAIN_SEED_LIST_H