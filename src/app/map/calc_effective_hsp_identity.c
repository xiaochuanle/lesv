#include "calc_effective_hsp_identity.h"

static int
find_next_trace_point(const char* qaln, const char* saln, const int aln_size, int* aln_idx)
{
    const int E = 20;
    int i = *aln_idx;
    while (i < aln_size) {
        while (i < aln_size && qaln[i] != saln[i]) ++i;
        if (i >= aln_size) break;
        int j = i + 1;
        while (j < aln_size && qaln[j] == saln[j]) ++j;
        if 
    }
}

static void
calc_trace_points(const char* qaln, const char* saln, const int aln_size, vec_int_pair* tp_list)
{
    static const int kSegLen = 20000;
    kv_clear(*tp_list);

}