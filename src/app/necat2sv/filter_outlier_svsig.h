#ifndef __FILTER_OUTLIER_SVSIG_H
#define __FILTER_OUTLIER_SVSIG_H

#include "sv_signature.h"

#ifdef __cplusplus
extern "C" {
#endif

void
find_outlier_svsig(SvSignature* svsig_array, int svsig_count);

void 
trim_svsig_by_cov_stats(SvSignature* svga, int svgc);

#ifdef __cplusplus
}
#endif

#endif // __FILTER_OUTLIER_SVSIG_H