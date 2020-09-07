#ifndef __CNS_OPTIONS_H
#define __CNS_OPTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const char* wrk_dir;
    const char* svr_group_dir;
    const char* db_dir;
    const char* db_title;

    int     num_threads;
    int     node_id;
    int     num_nodes;
    int     batch_size;
    int     use_batch_mode;
    double  ovlp_cov_perc;
    int     ovlp_cov_res;
    double  cns1_perc_identity;
    double  cns2_perc_identity;
    int     refine_ovlp;
    int     min_cov;
    int     min_size;

    int     memsc_kmer_size;
    int     memsc_kmer_window;
    int     memsc_mem_size;
    int     memsc_score;
} HbnProgramOptions;

#ifdef __cplusplus
}
#endif

#endif // __CNS_OPTIONS_H