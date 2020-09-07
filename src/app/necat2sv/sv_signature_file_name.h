#ifndef __SV_SIGNATURE_FILE_NAME_H
#define __SV_SIGNATURE_FILE_NAME_H

#ifdef __cplusplus
extern "C" {
#endif

void create_sv_signature_dir(const char* wrk_dir);

int load_subject_sv_signature_file_count(const char* wrk_dir);

void dump_subject_sv_signature_file_count(const char* wrk_dir, const int cnt);

char* make_subject_sv_signature_path(const char* wrk_dir, const int sid, char path[]);

int subject_sv_signature_is_created(const char* wrk_dir, const int sid);

void subject_sv_signature_make_created(const char* wrk_dir, const int sid);

#ifdef __cplusplus
}
#endif

#endif // __SV_SIGNATURE_FILE_NAME_H