#ifndef __MAP_RESULTS_H
#define __MAP_RESULTS_H

#include <stdio.h>

#include "../../corelib/hbn_aux.h"
#include "../../corelib/m4_record.h"
#include "../../corelib/seqdb.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const char* kSamVersion;

void
make_sam_rg_info(const int subject_id, kstring_t* rg);

#define DEFAULT_SAMPLE "mmsv"

void
print_sam_prolog(FILE* out, 
    const char* sam_version, 
    const char* prog_version, 
    const char* rg_info,
    const char* rg_sample,
    const CSeqDB* db,
    const int subject_id,
    int argc, 
    char* argv[]);

void
print_one_sam_result(
    const char* qname,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const char* sname,
    const int soff,
    const int send,
    const int ssize,
    const char* qaln,
    const char* saln,
    const int aln_size,
    const int map_score,
    const double perc_identity,
    const double eff_perc_identity,
    const BOOL dump_md,
    const char* rg_sample,
    kstring_t* out);

void
print_one_paf_result(
    const char* qname,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const char* sname,
    const int soff,
    const int send,
    const int ssize,
    const char* qaln,
    const char* saln,
    const int aln_size,
    const int map_score,
    const BOOL dump_cigar,
    const BOOL dump_md,
    kstring_t* out);

void
print_one_m4_result(
    const char* qname,
    const int qid,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const char* sname,
    const int sid,
    const int soff,
    const int send,
    const int ssize,
    const int map_score,
    const double perc_identity,
    const BOOL binary,
    kstring_t* line,
    kstring_t* out);

#ifdef __cplusplus
}
#endif

#endif // __MAP_RESULTS_H