#ifndef __HBN_WORD_FINDER_H
#define __HBN_WORD_FINDER_H

#include "../corelib/gapped_candidate.h"
#include "../corelib/seqdb.h"
#include "chain_dp.h"
#include "lookup_table.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	u64 hash;
    int context;
    int offset;
} DDFKmer;

typedef kvec_t(DDFKmer) vec_ddfk;

typedef struct {
    int context;
    int qoff;
    idx soff;
} DDFKmerMatch;

typedef kvec_t(DDFKmerMatch) vec_ddfkm;

typedef struct {
    DDFKmerMatch* km_array;
    int km_count;
} QueryKmerMatchInfo;

typedef struct {
    LookupTable* lktbl;
    const text_t* subject_vol;
    const text_t* query_vol;
    int qidx_from;
    int qidx_to;
    int thread_id;
    int num_threads;
    vec_ddfkm l_ddfkm_list;
    QueryKmerMatchInfo* g_kmi_array;
    int kmer_size;
    int kmer_window;
    int max_kmer_occ;
    const idx* global_soff_max_array;
} WordFinderThreadData;

WordFinderThreadData*
WordFinderThreadDataNew(
    LookupTable* lktbl,
    const text_t* subject_vol,
    const text_t* query_vol,
    int qidx_from,
    int qidx_to,
    int thread_id,
    int num_threads,
    QueryKmerMatchInfo* kmi_array,
    int kmer_size,
    int kmer_window,
    int max_kmer_occ);

WordFinderThreadData*
WordFinderThreadDataFree(WordFinderThreadData* data);

typedef struct {
    LookupTable* lktbl;
    const text_t* subject_vol;
    const text_t* query_vol;
    int num_threads;
    WordFinderThreadData** wftd_array;
    QueryKmerMatchInfo* kmi_array;
    int kmer_size;
    int kmer_window;
    int max_kmer_occ;
} HbnWordFindData;

HbnWordFindData*
HbnWordFindDataNew(LookupTable* lktbl,
    const text_t* subject_vol,
    const text_t* query_vol,
    int num_threads,
    int kmer_size,
    int kmer_window,
    int max_kmer_occ);

HbnWordFindData*
HbnWordFindDataFree(HbnWordFindData* data);

void
hbn_find_kmer_matches_mt(HbnWordFindData* data,
    const int qidx_from,
    const int qidx_to,
    const idx* soff_max_array);

#ifdef __cplusplus
}
#endif

#endif // __HBN_WORD_FINDER_H
