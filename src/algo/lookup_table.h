#ifndef __LOOKUP_TABLE_H
#define __LOOKUP_TABLE_H

#include "../corelib/seqdb.h"
#include "../ncbi_blast/setup/blast_sequence_blk.h"
#include "../ncbi_blast/setup/blast_query_info.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	u64 hash;
	u64 offset_list_idx;
	int cnt;
} LktblKmerInfo;

typedef struct {
	u64*	offset_list;
	u64 	offset_list_size;
	LktblKmerInfo* 	lki_array;
	u64 			lki_count;
} LookupTable;

LookupTable*
LookupTableFree(LookupTable* lktbl);

LookupTable*
build_lookup_table(const text_t* db,
    const int kmer_size,
    const int window_size,
    const int max_kmer_occ,
    const int num_threads);

#ifdef __cplusplus
}
#endif

#endif // __LOOKUP_TABLE_H