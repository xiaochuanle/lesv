#ifndef __SV_READ_FILE_NAME_H
#define __SV_READ_FILE_NAME_H

#ifdef __cplusplus
extern "C" {
#endif

void create_sv_read_dir(const char* wrk_dir);

char* make_query_vol_sv_read_path(const char* wrk_dir, const int vid, char path[]);

int query_vol_sv_read_is_created(const char* wrk_dir, const int vid);

void query_vol_sv_read_make_created(const char* wrk_dir, const int vid);

void dump_subject_sv_read_file_count(const char* wrk_dir, const int cnt);

int load_subject_sv_read_file_count(const char* wrk_dir);

char* make_subject_sv_read_path(const char* wrk_dir, const int sid, char path[]);

#ifdef __cplusplus
}
#endif

#endif // __SV_READ_FILE_NAME_H