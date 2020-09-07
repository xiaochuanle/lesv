#include "sv_read_group_file_name.h"

#include "../../corelib/hbn_aux.h"

#include <fstream>
#include <string>
#include <sstream>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>

using namespace std;

extern "C" 
void create_sv_read_group_dir(const char* wrk_dir)
{
    const char* path = wrk_dir;
    if ((access(path, F_OK) != 0)
        &&
        (mkdir(path, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", path, strerror(errno));
    }
}

extern "C"
char* make_subject_sv_read_group_path(const char* wrk_dir, const int sid, char path[])
{
    string sv_read_group_dir = wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    sprintf(path, "%s/subject_%s_sv_read_group", sv_read_group_dir.c_str(), sid_str);
    return path;
}

extern "C"
int subject_sv_read_group_is_created(const char* wrk_dir, const int sid)
{
    ostringstream os;
    os << wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    os << "/subject_" << sid_str << "_sv_read_group.created";
    string path = os.str();
    return access(path.c_str(), F_OK) == 0;
}

extern "C"
void subject_sv_read_group_make_created(const char* wrk_dir, const int sid)
{
    ostringstream os;
    os << wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    os << "/subject_" << sid_str << "_sv_read_group.created";
    string path = os.str();
    
    hbn_dfopen(out, path.c_str(), "w");
    hbn_fclose(out);
}

extern "C"
int subject_sv_read_group_is_corrected(const char* wrk_dir, const int sid)
{
    ostringstream os;
    os << wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    os << "/subject_" << sid_str << "_sv_read_group.corrected";
    string path = os.str();
    return access(path.c_str(), F_OK) == 0;    
}

extern "C"
void subject_sv_read_group_make_corrected(const char* wrk_dir, const int sid)
{
    ostringstream os;
    os << wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    os << "/subject_" << sid_str << "_sv_read_group.corrected";
    string path = os.str();
    
    hbn_dfopen(out, path.c_str(), "w");
    hbn_fclose(out);
}

extern "C"
void make_subject_corrected_normal_sv_read_path(const char* wrk_dir, const int sid, char path[])
{
    ostringstream os;
    os << wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    os << "/subject_" << sid_str << ".normal.cns.fasta";
    string cpp_path = os.str();
    memcpy(path, cpp_path.c_str(), cpp_path.size());
    path[cpp_path.size()] = '\0';    
}

extern "C"
void make_subject_sv_read_sam_path(const char* wrk_dir, const int sid, char path[])
{
    ostringstream os;
    os << wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    os << "/subject_" << sid_str << ".sam";
    string cpp_path = os.str();
    memcpy(path, cpp_path.c_str(), cpp_path.size());
    path[cpp_path.size()] = '\0';        
}

extern "C"
char* make_subject_outlier_fasta_path(const char* wrk_dir, const int sid, char path[])
{
    string sv_read_dir = wrk_dir;
    char sid_str[64];
    u64_to_fixed_width_string_r(sid, sid_str, HBN_DIGIT_WIDTH);
    sprintf(path, "%s/subject_%s.outlier.fasta", sv_read_dir.c_str(), sid_str);
    return path;
}

extern "C"
char* make_sv_read_header(const char* old_name, int qdir, int sid, int gid, int sfrom, int sto, kstring_t* new_name)
{
    ks_clear(*new_name);
    ksprintf(new_name, "%s_svr:%d:%d:%d:%d:%d", old_name, qdir, sid, gid, sfrom, sto);
    return ks_s(*new_name);
}

extern "C"
void
extract_info_from_sv_read_header(const char* name, int* qdir, int* sid, int* gid, int* sfrom, int* sto)
{
    const int n = strlen(name);
    int x = n;
    while (x > 2) {
        --x;
        if (name[x] =='r') {
            if (name[x-1] == 'v' && name[x-2] == 's') {
                x -= 2;
                break;
            }
        }
    }
    hbn_assert(x > 0);
    hbn_assert(name[x+3] == ':');
    x += 3;
    hbn_assert(name[x] == ':');
    ++x;
    if (qdir) *qdir = atoi(name + x);

    while (x < n) {
        if (name[x] == ':') break;
        ++x;
    }
    hbn_assert(x < n);
    ++x;
    if (sid) *sid = atoi(name + x);

    while (x < n) {
        if (name[x] == ':') break;
        ++x;
    }
    hbn_assert(x < n);
    ++x;
    if (gid) *gid = atoi(name + x);

    while (x < n) {
        if (name[x] == ':') break;
        ++x;
    }
    hbn_assert(x < n);
    ++x;
    if (sfrom) *sfrom = atoi(name + x);

    while (x < n) {
        if (name[x] == ':') break;
        ++x;
    }
    hbn_assert(x < n);
    ++x;
    if (sto) *sto = atoi(name + x);
}