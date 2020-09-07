#ifndef __REFINE_ALIGN_H
#define __REFINE_ALIGN_H

#include "chain_dp.h"
#include "dalign.h"
#include "edlib_wrapper.h"
#include "ksw2.h"
#include "ksw2_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

int
refine_align(Ksw2Data* ksw,
    const u8* query,
    const int query_align_from,
    const int query_align_to,
    const int query_length,
    const u8* subject,
    const int subject_align_from,
    const int subject_align_to,
    const int subject_length,
    const char* org_qaln,
    const char* org_saln,
    const int org_align_size,
    kstring_t* qaln,
    kstring_t* saln);

#ifdef __cplusplus
}
#endif

#endif // __REFINE_ALIGN_H