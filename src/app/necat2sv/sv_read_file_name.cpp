#include "sv_read_file_name.h"

#include "../../corelib/hbn_aux.h"

#include <fstream>
#include <string>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>

using namespace std;

extern "C" 
void create_sv_read_dir(const char* wrk_dir)
{
    const char* path = wrk_dir;
    if ((access(path, F_OK) != 0)
        &&
        (mkdir(path, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", path, strerror(errno));
    }
}

extern "C"
char* make_query_vol_sv_read_path(const char* wrk_dir, const int vid, char path[])
{
    char vid_str[64];
    u64_to_fixed_width_string_r(vid, vid_str, HBN_DIGIT_WIDTH);
    sprintf(path, "%s/Q%s.sv_read", wrk_dir, vid_str);
    return path;
}

extern "C"
int query_vol_sv_read_is_created(const char* wrk_dir, const int vid)
{
    char path[HBN_MAX_PATH_LEN];
    make_query_vol_sv_read_path(wrk_dir, vid, path);
    strcat(path, ".created");
    return access(path, F_OK) == 0;
}

extern "C"
void query_vol_sv_read_make_created(const char* wrk_dir, const int vid)
{
    char path[HBN_MAX_PATH_LEN];
    make_query_vol_sv_read_path(wrk_dir, vid, path);
    strcat(path, ".created");
    hbn_dfopen(out, path, "w");
    hbn_fclose(out);
}

extern "C"
void dump_subject_sv_read_file_count(const char* wrk_dir, const int cnt)
{
    string path = wrk_dir;
    path += "/subject_sv_read_file_count";
    hbn_dfopen(out, path.c_str(), "w");
    fprintf(out, "%d\n", cnt);
    hbn_fclose(out);
}

extern "C"
int load_subject_sv_read_file_count(const char* wrk_dir)
{
    string path = wrk_dir;
    path += "/subject_sv_read_file_count";
    int cnt;
    ifstream in(path);
    in >> cnt;
    in.close();
    return cnt;
}

extern "C"
char* make_subject_sv_read_path(const char* wrk_dir, const int sid, char path[])
{
    string sv_read_dir = wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    sprintf(path, "%s/subject_%s.sv_read", sv_read_dir.c_str(), sid_str);
    return path;
}