#ifndef __CNS_READ_HEADER_H
#define __CNS_READ_HEADER_H

#include "kstring.h"
#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

void
make_cns_read_header(const char* old_name,
    int num_extended_sr,
    int num_added_aln,
    int cns_subseq_offset,
    int cns_subseq_end,
    int cns_seq_size,
    int raw_seq_size,
    kstring_t* new_name);

void
extract_info_from_cns_read_header(const char* name,
    int* num_extended_sr,
    int* num_added_aln,
    int* cns_subseq_offset,
    int* cns_subseq_end,
    int* cns_seq_size,
    int* raw_seq_size);

#ifdef __cplusplus
}
#endif

#endif // __CNS_READ_HEADER_H