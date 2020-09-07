#include "cns_read_header.h"

void
make_cns_read_header(const char* old_name,
    int num_extended_sr,
    int num_added_aln,
    int cns_subseq_offset,
    int cns_subseq_end,
    int cns_seq_size,
    int raw_seq_size,
    kstring_t* new_name)
{
    ks_clear(*new_name);
    ksprintf(new_name, "%s_cns:%d:%d:%d:%d:%d:%d",
        old_name,
        num_extended_sr,
        num_added_aln,
        cns_subseq_offset,
        cns_subseq_end,
        cns_seq_size,
        raw_seq_size);
}

void
extract_info_from_cns_read_header(const char* name,
    int* num_extended_sr,
    int* num_added_aln,
    int* cns_subseq_offset,
    int* cns_subseq_end,
    int* cns_seq_size,
    int* raw_seq_size)
{
    int n = strlen(name);
    int x = n;
    while (x > 2) {
        --x;
        if (name[x] == 's') {
            if (name[x-1] == 'n' && name[x-2] == 'c') {
                x -= 2;
                break;
            }
        }
    }
    hbn_assert(x > 0);
    hbn_assert(name[x+3] == ':');

    int t;
    x += 4;
    hbn_assert(x < n);
    t = atoi(name + x);
    if (num_extended_sr) *num_extended_sr = t;

    while (x < n && name[x] != ':') ++x;
    hbn_assert(name[x] == ':');
    ++x;
    hbn_assert(x < n);
    t = atoi(name + x);
    if (num_added_aln) *num_added_aln = t;

    while (x < n && name[x] != ':') ++x;
    hbn_assert(name[x] == ':');
    ++x;
    hbn_assert(x < n);
    t = atoi(name + x);
    if (cns_subseq_offset) *cns_subseq_offset = t;   

    while (x < n && name[x] != ':') ++x;
    hbn_assert(name[x] == ':');
    ++x;
    hbn_assert(x < n);
    t = atoi(name + x);
    if (cns_subseq_end) *cns_subseq_end = t;   

    while (x < n && name[x] != ':') ++x;
    hbn_assert(name[x] == ':');
    ++x;
    hbn_assert(x < n);
    t = atoi(name + x);
    if (cns_seq_size) *cns_seq_size = t;   

    while (x < n && name[x] != ':') ++x;
    hbn_assert(name[x] == ':');
    ++x;
    hbn_assert(x < n);
    t = atoi(name + x);
    if (raw_seq_size) *raw_seq_size = t;   
}