#ifndef __RESOVLE_BROKEN_MAP_H
#define __RESOVLE_BROKEN_MAP_H

#include "map_aux.h"
#include "../../corelib/m4_record.h"

#ifdef __cplusplus
extern "C" {
#endif

int
resolve_broken_map(HbnTracebackData* tbck_data,
    BlastHSP* hsp_array,
    int hsp_cnt,
    const char* read_name,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_length,
    const int subject_id,
    const char* subject_name,
    const u8* subject,
    const int subject_length,
    const double min_perc_identity,
    int* qdir_,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    int* score_,
    double* perc_identity_,
    double* eff_perc_identity_,
    kstring_t* qaln,
    kstring_t* saln);

#ifdef __cplusplus
}
#endif

#endif // __RESOVLE_BROKEN_MAP_H