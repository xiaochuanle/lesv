#ifndef __MAP_ONE_READ_H
#define __MAP_ONE_READ_H

#include "map_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

void
asm_map_one_read(AsmMapThreadData* data, int sr_info_idx);

#ifdef __cplusplus
}
#endif

#endif // __MAP_ONE_READ_H