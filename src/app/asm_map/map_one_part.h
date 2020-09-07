#ifndef __MAP_ONE_PART_H
#define __MAP_ONE_PART_H

#include "map_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

void
asm_map_one_part(const char* psr_dir,
    const HbnProgramOptions* opts, 
    RawReadReader* raw_reads, 
    const int pid,
    const int batch_size,
    const BOOL use_batch_mode);

#ifdef __cplusplus
}
#endif

#endif  // __MAP_ONE_PART_H