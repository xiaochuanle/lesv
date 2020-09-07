#ifndef __APPROX_SEMI_GAPPED_ALIGN
#define __APPROX_SEMI_GAPPED_ALIGN

#include "ksw2_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int qid;
    int sid;
    Ksw2Data* ksw;
    const u8* fwd_query;
    int query_length;
    const u8* fwd_subject;
    int subject_length;
    vec_u8 rev_query;
    vec_u8 rev_subject;
    kstring_t qaux;
    kstring_t saux;
    kstring_t fwd_qaln;
    kstring_t fwd_saln;
    kstring_t rev_qaln;
    kstring_t rev_saln;
    kstring_t bp_qaln;
    kstring_t bp_saln;
    int f_qoff, f_qend;
    int f_soff, f_send;
    int r_qoff, r_qend;
    int r_soff, r_send;
    double perc_identity, eff_perc_identity;
    kstring_t qaln;
    kstring_t saln;
    int qoff, qend;
    int soff, send;
} ApproxSemeGappedAlignData;

#define ApproxSemeGappedAlignData_DumpAlignPos(output_func_, stream_, asgad_) \
    output_func_(stream_, "fwd: [%d, %d, %d] x [%d, %d, %d]\n" \
        "rev: [%d, %d, %d] x [%d, %d, %d]\n" \
        "frev: [%d, %d, %d] x [%d, %d, %d]\n", \
        (asgad_)->f_qoff, \
        (asgad_)->f_qend, \
        (asgad_)->query_length, \
        (asgad_)->f_soff, \
        (asgad_)->f_send, \
        (asgad_)->subject_length, \
        (asgad_)->r_qoff, \
        (asgad_)->r_qend, \
        (asgad_)->query_length, \
        (asgad_)->r_soff, \
        (asgad_)->r_send, \
        (asgad_)->subject_length, \
        (asgad_)->query_length - (asgad_)->r_qend, \
        (asgad_)->query_length - (asgad_)->r_qoff, \
        (asgad_)->query_length, \
        (asgad_)->subject_length - (asgad_)->r_send, \
        (asgad_)->subject_length - (asgad_)->r_soff, \
        (asgad_)->subject_length)

#define ApproxSemeGappedAlignData_Dump(output_func_, stream_, asgad_) \
    output_func_(stream_, "[%d, %d, %d] x [%d, %d, %d], perc_identity = %g, eff_perc_identity = %g\n", \
        (asgad_)->qoff, \
        (asgad_)->qend, \
        (asgad_)->query_length, \
        (asgad_)->soff, \
        (asgad_)->send, \
        (asgad_)->subject_length, \
        (asgad_)->perc_identity, \
        (asgad_)->eff_perc_identity)


ApproxSemeGappedAlignData*
ApproxSemeGappedAlignDataNew();

ApproxSemeGappedAlignData*
ApproxSemeGappedAlignDataFree(ApproxSemeGappedAlignData* data);

void
ApproxSemeGappedAlignData_Init(ApproxSemeGappedAlignData* data,
    const u8* query, const int query_length,
    const u8* subject, const int subject_length);

BOOL 
ApproxSemeGappedAlignData_SemiGappedAlign(ApproxSemeGappedAlignData* data,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length);

void
ApproxSemeGappedAlignData_UnitTest(const char* fasta_path);

#ifdef __cplusplus
}
#endif

#endif // __APPROX_SEMI_GAPPED_ALIGN
