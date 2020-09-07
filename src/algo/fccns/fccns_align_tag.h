#ifndef __FCCNS_ALIGN_TAG_H
#define __FCCNS_ALIGN_TAG_H

#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DEFAULT_CNS_WEIGHT  1.0

typedef u16 align_tag_delta_t;
#define ALIGN_TAG_MAX_DELTA U16_MAX

typedef struct {
    double weight;
    int t_pos;
    int p_t_pos;
    align_tag_delta_t delta;
    align_tag_delta_t p_delta;
    char q_base;
    char p_q_base;
} AlignTag;

#define dump_align_tag(output_func, stream, tag) \
    output_func(stream, "weight = %g, t_pos = %d, p_t_pos = %d, delta = %d, p_delta = %d, q_base = %c, p_q_base = %c\n", \
        (tag).weight, \
        (tag).t_pos, \
        (tag).p_t_pos, \
        (tag).delta, \
        (tag).p_delta, \
        (tag).q_base, \
        (tag).p_q_base)

#define align_tag_plink_eq(a, b) \
    ((a).p_t_pos == (b).p_t_pos && (a).p_delta == (b).p_delta && (a).p_q_base == (b).p_q_base)

typedef kvec_t(AlignTag) vec_align_tag;

void ks_introsort_align_tag_lt(size_t n, AlignTag* a);

void
make_align_tags_from_ovlp(const char* qaln,
    const char* taln,
    const size_t aln_size,
    const int qoff,
    const int qend,
    const int toff,
    const int tend,
    const double weight,
    vec_align_tag* align_tag_list);

#ifdef __cplusplus
}
#endif

#endif // __FCCNS_ALIGN_TAG_H