#ifndef __SV_SIGNATURE_H
#define __SV_SIGNATURE_H

#include "../../corelib/kvec.h"
#include "../../ncbi_blast/setup/gapinfo.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int qid;
    int qdir;
    int qfrom, qto;
    int fqfrom, fqto;
    int qsize;
    int sfrom, sto;
    int fsfrom, fsto;
    int ssize;
    EGapAlignOpType type;
} SvSignature;

typedef kvec_t(SvSignature) vec_sv_sig;

void sort_svsig_sfrom_lt(size_t n, SvSignature* a);

#define dump_svsig(output_func, stream, svpos, qname_) do { \
    output_func(stream, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", \
        (svpos).qid, \
        (svpos).qdir, \
        (svpos).qfrom, \
        (svpos).qto, \
        (svpos).fqfrom, \
        (svpos).fqto, \
        (svpos).qsize, \
        (svpos).sfrom, \
        (svpos).sto, \
        (svpos).fsfrom, \
        (svpos).fsto, \
        (svpos).ssize, \
        (svpos).type); \
    const char* __qname = (qname_); \
    if (__qname) output_func(stream, "\t%s", __qname); \
    output_func(stream, "\n"); \
} while(0)

void
sread_sv_signature(const char* in, SvSignature* svsig);

SvSignature*
load_sv_signature_array(const char* path, int* sv_signature_count, const EGapAlignOpType type);

#ifdef __cplusplus
}
#endif

#endif // __SV_SIGNATURE_H