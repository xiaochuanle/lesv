#ifndef __SV_READ_GROUP_FILE_NAME_H
#define __SV_READ_GROUP_FILE_NAME_H

#include "../../corelib/kstring.h"

#ifdef __cplusplus
extern "C" {
#endif

void create_sv_read_group_dir(const char* wrk_dir);

char* make_subject_sv_read_group_path(const char* wrk_dir, const int sid, char path[]);

int subject_sv_read_group_is_created(const char* wrk_dir, const int sid);

void subject_sv_read_group_make_created(const char* wrk_dir, const int sid);

void subject_sv_read_group_make_corrected(const char* wrk_dir, const int sid);

int subject_sv_read_group_is_corrected(const char* wrk_dir, const int sid);

void make_subject_corrected_normal_sv_read_path(const char* wrk_dir, const int sid, char path[]);

void make_subject_sv_read_sam_path(const char* wrk_dir, const int sid, char path[]);

char* make_subject_outlier_fasta_path(const char* wrk_dir, const int sid, char path[]);

char* make_sv_read_header(const char* old_name, int qdir, int sid, int gid, int sfrom, int sto, kstring_t* new_name);

void extract_info_from_sv_read_header(const char* name, int* qdir, int* sid, int*gid, int* sfrom, int* sto);

#ifdef __cplusplus
}
#endif

#endif // __SV_READ_GROUP_FILE_NAME_H