#ifndef __CNS_ONE_GROUP_H
#define __CNS_ONE_GROUP_H

#include "cns_aux.h"
#include "sv_read_group.h"

#ifdef __cplusplus
extern "C" {
#endif

int 
cns_next_group(CnsThreadData* data, 
    FILE* sv_read_group_file,
    const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    FILE* packed_read_file,
    const int subject_id);

int 
cns_one_group_mt(CnsThreadData** data, 
    FILE* sv_read_group_file,
    const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    FILE* packed_read_file,
    const int subject_id);

void
map_sv_read_group(CnsThreadData* data, SvReadGroup* group);

void
cns_suject_svr_groups_st(int argc, char* argv[],
    HbnProgramOptions* opts,
    CSeqDB* db,
    CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    FILE* pac_file,
    const int sid);

#ifdef __cplusplus
}
#endif

#endif // __CNS_ONE_GROUP_H