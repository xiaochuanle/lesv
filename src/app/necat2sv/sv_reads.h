#ifndef __SV_READS_H
#define __SV_READS_H

#include <stdlib.h>
#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int query_id;
    int subject_id;
    int qdir, qoff, qend, qsize;
    int soff, send;
    int dist;
    int dual_qdir, dual_qoff, dual_qend;
    int dual_soff, dual_send;
    int dual_dist;
} SvRead;

typedef kvec_t(SvRead) vec_sv_read;

#define init_sv_read(read_) ((read_).query_id = (read_).subject_id = (read_).dual_qdir = -1)

#define dump_sv_read(output_func, stream, read, read_name_) do {\
    output_func(stream, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", \
        (read).query_id, \
        (read).qdir, \
        (read).qoff, \
        (read).qend, \
        (read).qsize, \
        (read).subject_id, \
        (read).soff, \
        (read).send, \
        (read).dist); \
    if ((read).dual_qdir != -1) { \
        output_func(stream, "\t%d\t%d\t%d\t%d\t%d\t%d", \
            (read).dual_qdir, \
            (read).dual_qoff, \
            (read).dual_qend, \
            (read).dual_soff, \
            (read).dual_send, \
            (read).dual_dist); \
    } \
    const char* __read_name = (read_name_); \
    if (__read_name) output_func(stream, "\t%s", __read_name); \
    output_func(stream, "\n"); \
} while(0)

void
sread_sv_read(const char* in, SvRead* rp);

SvRead*
load_sv_read_array(const char* path, int* sv_read_count);
 
#ifdef __cplusplus
}
#endif

#endif // __SV_READS_H