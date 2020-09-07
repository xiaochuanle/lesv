#ifndef __CHAIN_AND_EXTEND_KMER_MATCHES_H
#define __CHAIN_AND_EXTEND_KMER_MATCHES_H

#include "../../algo/hbn_word_finder.h"
#include "hbn_options.h"

#ifdef __cplusplus
extern "C" {
#endif

void
chain_and_extend_km_mt(HbnWordFindData* word_data,
    const HbnProgramOptions* opts,
    CSeqDB* query_vol,
    CSeqDB* subject_vol,
    const int qidx_from,
    const int qidx_to,
    FILE* out,
    FILE* backup_out);

#ifdef __cplusplus
}
#endif

#endif // __CHAIN_AND_EXTEND_KMER_MATCHES_H