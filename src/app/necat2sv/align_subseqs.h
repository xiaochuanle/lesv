#ifndef __ALIGN_SUBSEQS_H
#define __ALIGN_SUBSEQS_H

#include "../../algo/hbn_traceback.h"
#include "../../corelib/raw_reads_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    kstring_t qaux;
    kstring_t saux;
    const char* qas;
    const char* qae;
    const char* sas;
    const char* sae;
    int qb;
    int qe;
    int sb;
    int se;
    int dist;
    int score;
    double perc_identity;
    double eff_perc_identity;
} SubseqAlignResult;

SubseqAlignResult*
SubseqAlignResultNew();

SubseqAlignResult*
SubseqAlignResultFree(SubseqAlignResult* result);

BOOL
align_and_refine_subseq_with_edlib(RawReadReader* query_reader,
    RawReadReader* subject_reader,
    HbnTracebackData* tbck_data,
    const int max_dist,
    int query_gid, 
    int query_dir,
    int query_from,
    int query_to,
    int subject_gid,
    int subject_from,
    int subject_to,
    SubseqAlignResult* result);

BOOL
align_and_refine_subseq_with_ksw(RawReadReader* query_reader,
    RawReadReader* subject_reader,
    HbnTracebackData* tbck_data,
    const int max_dist,
    int query_gid, 
    int query_dir,
    int query_from,
    int query_to,
    int subject_gid,
    int subject_from,
    int subject_to,
    SubseqAlignResult* result);

#ifdef __cplusplus
}
#endif

#endif // __ALIGN_SUBSEQS_H