#ifndef __SORT_SR_HIT_SEEDS_H
#define __SORT_SR_HIT_SEEDS_H

#include "hbn_word_finder.h"

#ifdef __cplusplus
extern "C" {
#endif

void
sort_ddfkm_context_lt(size_t n, DDFKmerMatch* a);

void
sort_ddfkm_soff_lt(size_t n, DDFKmerMatch* a);

void
sort_ddfk_hash_lt(size_t n, DDFKmer* a);

void
sort_chain_seed_soff_lt(size_t n, ChainSeed* a);

#ifdef __cplusplus
}
#endif

#endif // __SORT_SR_HIT_SEEDS_H