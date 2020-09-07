#include "partition_mt.h"

#include <pthread.h>

#include "cstr_util.h"

void
make_partition_name(const char* data_dir, const char* prefix, const int pid, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (prefix) strcat(path, prefix);
    char buf[64];
    u64_to_fixed_width_string_r(pid, buf, HBN_DIGIT_WIDTH);
    strcat(path, buf);
}

void
dump_partition_count(const char* data_dir, const char* prefix, const int np)
{
    char path[HBN_MAX_PATH_LEN];
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (prefix) {
        strcat(path, prefix);
        strcat(path, ".");
    }
    strcat(path, "np");
    hbn_dfopen(out, path, "w");
    fprintf(out, "%d\n", np);
    hbn_fclose(out);
}

int
load_partition_count(const char* data_dir, const char* prefix)
{
    char path[HBN_MAX_PATH_LEN];
    if (prefix) {
        sprintf(path, "%s/%s.np", data_dir, prefix);
    } else {
        sprintf(path, "%s/np", data_dir);
    }
    int np;
    hbn_dfopen(in, path, "r");
    HBN_SCANF(fscanf, in, 1, "%d", &np);
    hbn_fclose(in);
    return np;    
}

void* load_part_records(const char* path, const size_t record_size, size_t* n_record)
{
    size_t s = hbn_file_size(path);
    if (!s) return NULL;
    hbn_assert(s % record_size == 0, "s = %zu, n = %zu", s, record_size);
    void* a = malloc(s);
    s /= record_size;
    hbn_dfopen(in, path, "rb");
    hbn_fread(a, record_size, s, in);
    hbn_fclose(in);
    *n_record = s;
    return a;    
}

/////////////////////////

static size_t 
load_next_batch_records(void* array, const size_t size, size_t count, FILE* in)
{
    size_t n = fread(array, size, count, in);
    return n;
}

#define id_in_range(id_, L_, R_) ((id_) >= (L_) && (id_) < (R_))

static void
process_one_batch(const char* part_wrk_dir,
    void* a,
    size_t n,
    const size_t record_size,
    const int batch_size,
    const int num_batches,
    const int num_threads,
    qid_extract_func get_qid,
    sid_extract_func get_sid,
    change_record_roles_func change_roles,
    normolize_sdir_func normalise_sdir,
    sort_record_mt_func sort_record_mt)
{
    if (n == 0) return;

    char ne[record_size];
    size_t m = n;
    if (change_roles) {
        for (size_t i = 0; i < n; ++i) {
            void* e = a + record_size * i;
	        //int qid = get_qid(e);
	        //int sid = get_sid(e);
            change_roles(e, ne);
            if (normalise_sdir) {
                normalise_sdir(e);
                normalise_sdir(ne);
            }
            memcpy(a + m * record_size, ne, record_size);
            ++m;
        }
    }
    HBN_LOG("n = %zu, m = %zu", n, m);
    n = m;
    (*sort_record_mt)(n, a, num_threads);
    for (size_t i = 0; i < n - 1; ++i) {
        void* e = a + i * record_size;
        void* e1 = e + record_size;
        int sid = get_sid(e);
        int sid1 = get_sid(e1);
        hbn_assert(sid <= sid1);
    }
    
    size_t i = 0;
    char path[HBN_MAX_PATH_LEN];
    while (i < n) {
        void* e = a + i * record_size;
        int sid = get_sid(e);
	    int qid = get_qid(e);
        int pid = sid / batch_size;
        hbn_assert(pid < num_batches, "i = %zu, record_size = %zu, qid = %d, sid = %d, pid = %d, batch_size = %d, num_batches = %d",
		i, record_size, qid, sid, pid, batch_size, num_batches);
        int max_sid = (pid + 1) * batch_size;
        size_t j = i + 1;
        while (j < n) {
            e += record_size;
            int sid1 = get_sid(e);
            if (sid1 >= max_sid) break;
            ++j;
        }
        make_partition_name(part_wrk_dir, DEFAULT_PART_PREFIX, pid, path);
        HBN_LOG("dump to file %s", path);
        hbn_dfopen(out, path, "ab");
        void* ap = a + i * record_size;
        m = j - i;
        HBN_LOG("dump %zu --- %zu", i, j);
        hbn_fwrite(ap, record_size, m, out);
        hbn_fclose(out);
        i = j;
    }
}

static void
create_partition_files(const char* part_wrk_dir, const int pid_from, const int pid_to)
{
    char path[HBN_MAX_PATH_LEN];
    for (int i = pid_from; i < pid_to; ++i) {
        make_partition_name(part_wrk_dir, DEFAULT_PART_PREFIX, i, path);
        hbn_dfopen(out, path, "wb");
        hbn_fclose(out);
    }
}

void
partition_record_mt(const char* part_wrk_dir,
    const char* record_path,
    const int num_batches,
    const int batch_size,
    const int num_threads,
    const size_t record_size,
    size_t work_space,
    qid_extract_func get_qid,
    sid_extract_func get_sid,
    change_record_roles_func change_roles,
    normolize_sdir_func normalise_sdir,
    sort_record_mt_func sort_record_mt)
{
    create_partition_files(part_wrk_dir, 0, num_batches);
    size_t total_size = hbn_file_size(record_path);
    if (work_space > total_size * 2) work_space = total_size * 2;
    hbn_assert((total_size % record_size) == 0);
    if (record_size == 0) {
        HBN_LOG("file [%s] is empty.", record_path);
        return;
    }
    char total_size_str[64];
    u64_to_string_datasize(total_size, total_size_str);
    size_t process_size = 0;
    char process_size_str[64];
    size_t curr_size = 0;
    char curr_size_str[64];

    size_t N2 = work_space / record_size;
    size_t N = N2 / 2;
    hbn_dfopen(in, record_path, "rb");
    void* a = malloc(work_space);
    while (1) {
        size_t n = load_next_batch_records(a, record_size, N, in);
        if (n == 0) break;
        curr_size = n * record_size;
        u64_to_string_datasize(curr_size, curr_size_str);
        HBN_LOG("process %s data", curr_size_str);
        process_one_batch(part_wrk_dir,
            a,
            n,
            record_size,
            batch_size,
            num_batches,
            num_threads,
            get_qid,
            get_sid,
            change_roles,
            normalise_sdir,
            sort_record_mt);
        process_size += curr_size;
        u64_to_string_datasize(process_size, process_size_str);
        HBN_LOG("%s/%s data processed", process_size_str, total_size_str);
    }
    free(a);
}
