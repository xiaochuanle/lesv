#ifndef __TRF_ARRAY_H
#define __TRF_ARRAY_H

#include "../../corelib/name2id_map.h"
#include "../../corelib/seqdb.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char* trf_array;
    CSeqDBInfo dbinfo;
    const CSeqInfo* seqinfo_array;
} TrfArray;

TrfArray*
TrfArrayBuild(const char* trf_path, const char* db_dir, const CSeqInfo* seqinfo_array);

TrfArray*
TrfArrayFree(TrfArray* array);

BOOL
fall_in_trf_subseq(TrfArray* array, int sid, size_t from, size_t to);

#ifdef __cplusplus
}
#endif

#endif // __TRF_ARRAY_H