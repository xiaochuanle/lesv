#ifndef __SV_READ_GROUP_H
#define __SV_READ_GROUP_H

#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/c_ncbi_blast_aux.h"
#include "../../ncbi_blast/setup/gapinfo.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int local_id;
    int global_id;
    const char* name;
    int size;
    size_t offset_in_sv_read;
    int raw_seq_size;
    int raw_seq_from;
    int raw_seq_to;
    int fsfrom;
    int fsto;
    int fsqdir;
} SvReadInfo;

typedef struct {
    int group_id;
    int subject_id;
    EGapAlignOpType sv_type;
    const CSeqInfo* raw_read_info_array;
    const char* raw_read_name_array;
    int sv_read_group_size;
    SvReadInfo* svri_array;
    vec_u8 fwd_sv_read;
    vec_u8 rev_sv_read;
} SvReadGroup;

SvReadGroup*
SvReadGroupLoad(const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    const int subject_id,
    FILE* sv_read_group_file,
    FILE* pac_file);

SvReadGroup*
SvReadGroupNew(const CSeqInfo* raw_read_info_array,
    const char* raw_read_name_array,
    const int group_id,
    const int subject_id,
    EGapAlignOpType sv_type,
    const int num_reads);

SvReadGroup*
SvReadGroupFree(SvReadGroup* group);

const u8* 
SvReadGroup_ExtractSeq(SvReadGroup* group, const int local_read_id, const int strand);

int
SvReadGroup_ReadSize(SvReadGroup* group, const int local_read_id);

int 
SvReadGroup_GlobalReadId(SvReadGroup* group, const int local_read_id);

const char*
SvReadGroup_ReadName(SvReadGroup* group, const int local_read_id);

#ifdef __cplusplus
}
#endif

#endif  // __SV_READ_GROUP_H