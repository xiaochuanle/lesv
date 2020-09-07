#ifndef __CORRECT_ONE_PART_H
#define __CORRECT_ONE_PART_H

#include "cns_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

void
correct_one_part(const HbnProgramOptions* opts, RawReadReader* raw_reads, const int pid);

#ifdef __cplusplus
}
#endif

#endif  // __CORRECT_ONE_PART_H