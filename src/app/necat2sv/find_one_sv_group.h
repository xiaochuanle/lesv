#ifndef __FIND_ONE_SV_GROUP_H
#define __FIND_ONE_SV_GROUP_H

#include "sv_signature.h"
#include "../../corelib/khash.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

BOOL 
find_next_ins_group(SvSignature* sga, 
    int sgc, 
    int* next_i, 
    int* group_id, 
    FILE* out,
    void* added_id_set);

BOOL
find_next_ins_group_relax(SvSignature* sga,
    int sgc,
    int* next_i,
    int* group_id,
    FILE* out,
    void* added_id_set);

BOOL 
find_next_del_group(SvSignature* sga, 
    int sgc, 
    int* next_i, 
    int* group_id, 
    FILE* out,
    void* added_id_set);

BOOL
find_next_del_group_relax(SvSignature* sga,
    int sgc,
    int* next_i,
    int* group_id,
    FILE* out,
    void* added_id_set);

#ifdef __cplusplus
}
#endif

#endif // __FIND_ONE_SV_GROUP_H