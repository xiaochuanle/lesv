#ifndef __HBN_EXTEND_SUBSEQ_HIT_DIFF_H
#define __HBN_EXTEND_SUBSEQ_HIT_DIFF_H

#include "hbn_extend_subseq_hit.h"

#ifdef __cplusplus
extern "C" {
#endif

void
hbn_extend_query_subseq_hit_list_diff(HbnSubseqHit* subseq_hit_array,
    int subseq_hit_count,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_id,
    const int query_length,
    const text_t* db,
    const HbnProgramOptions* opts,
    HbnSubseqHitExtnData* data,
    BlastHitList* hit_list,
    HbnHSPResults* results);

#ifdef __cplusplus
}
#endif

#endif // __HBN_EXTEND_SUBSEQ_HIT_DIFF_H