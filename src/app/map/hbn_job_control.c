#include "hbn_job_control.h"

const char*
make_qi_vs_sj_results_path(const char* wrk_dir, const char* stage, const int qi, const int sj, char path[])
{
    path[0] = '\0';
    char qi_buf[64];
    char sj_buf[64];
    u64_to_fixed_width_string_r(qi, qi_buf, HBN_DIGIT_WIDTH);
    u64_to_fixed_width_string_r(sj, sj_buf, HBN_DIGIT_WIDTH);

    if (wrk_dir) {
        sprintf(path, "%s/%s/stage%s_Q%s_D%s", wrk_dir, stage, stage, qi_buf, sj_buf);
    } else {
        sprintf(path, "%s/stage%s_Q%s_D%s", stage, stage, qi_buf, sj_buf);
    }

    return path;
}

FILE*
open_qi_vs_sj_results_file(const char* wrk_dir, const char* stage, const int qi, const int sj)
{
    char path[HBN_MAX_PATH_LEN];
    make_qi_vs_sj_results_path(wrk_dir, stage, qi, sj, path);
    hbn_dfopen(out, path, "w");
    return out;
}

BOOL 
qi_vs_sj_is_mapped(const char* wrk_dir, const char* stage, const int qi, const int sj)
{
    char path[HBN_MAX_PATH_LEN];
    make_qi_vs_sj_results_path(wrk_dir, stage, qi, sj, path);
    strcat(path, ".mapped");
    return (access(path, F_OK) == 0);
}

void
qi_vs_sj_make_mapped(const char* wrk_dir, const char* stage, const int qi, const int sj)
{
    char path[HBN_MAX_PATH_LEN];
    make_qi_vs_sj_results_path(wrk_dir, stage, qi, sj, path);
    strcat(path, ".mapped");
    hbn_dfopen(out, path, "w");
    hbn_fclose(out);    
}

BOOL
all_vs_sj_is_mapped(const char* wrk_dir, 
    const char* stage, 
    const int qi_start,
    const int num_query_vols,
    const int sj,
    const int node_id,
    const int num_nodes)
{
    for (int i = qi_start + node_id; i < num_query_vols; i += num_nodes) {
        if (!qi_vs_sj_is_mapped(wrk_dir, stage, i, sj)) return FALSE;
    }
    return TRUE;
}

void
merge_qi_vs_sj_results(const char* wrk_dir, const char* stage, const int qi, const int sj, FILE* out)
{
    char path[HBN_MAX_PATH_LEN];
    make_qi_vs_sj_results_path(wrk_dir, stage, qi, sj, path);
    hbn_dfopen(qi_vs_sj_out, path, "r");
    const int kBufLen = 2028;
    char buf[kBufLen];
    while (1) {
        int n = fread(buf, 1, kBufLen, qi_vs_sj_out);
        fwrite(buf, 1, n, out);
        if (n < kBufLen) break;
    }
    hbn_fclose(qi_vs_sj_out);
}

void
merge_all_vs_sj_results(const char* wrk_dir, 
    const char* stage, 
    const int qi_start,
    const int num_query_vols,
    const int sj,
    const int node_id,
    const int num_nodes,
    FILE* out)
{
    for (int i = qi_start + node_id; i < num_query_vols; i += num_nodes) {
        merge_qi_vs_sj_results(wrk_dir, stage, i, sj, out);
    }
}