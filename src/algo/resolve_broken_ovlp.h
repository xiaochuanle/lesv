#ifndef __RESOLVE_BROKEN_OVLP_H
#define __RESOLVE_BROKEN_OVLP_H

#include "cns_extend_chain_seed_list.h"
#include "../corelib/m4_record.h"

#ifdef __cplusplus
extern "C" {
#endif

BOOL 
is_correct_ovlp_positions(const int qb, const int qe, const int qs,
    const int sb, const int se, const int ss);

BOOL
is_overhang_ovlp(const int qb, const int qe, const int qs, 
    const int sb, const int se, const int ss);

int
resolve_broken_ovlp(HbnTracebackData* tbck_data, 
    SupportReadHitFindData* hit_finder, 
    vec_chain_seed* chain_seed_list,
    HbnInitHit* hit_array, 
    int hit_count,
    const char* read_name,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_length,
    const int subject_id,
    const char* subject_name,
    const u8* subject,
    const int subject_length,
    const int min_cov_res,
    const double min_perc_coverage,
    const double min_perc_identity,
    int* qdir_,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* perc_identity_,
    double* eff_perc_identity_,
    kstring_t* qaln,
    kstring_t* saln);

#ifdef __cplusplus
}
#endif

#endif // __RESOLVE_BROKEN_OVLP_H