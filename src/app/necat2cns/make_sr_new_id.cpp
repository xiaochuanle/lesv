#include "../../corelib/seqdb.h"
#include "../../corelib/partition_mt.h"
#include "../../corelib/gapped_candidate.h"

#include <utility>

using namespace std;

static void
s_add_one_part_seq_ids(const char* cns_sr_dir, int pid, int* local_id, pair<int, int>* id_map_array)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(cns_sr_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".cns_seq_id");
    size_t fs = hbn_file_size(path);
    hbn_assert((fs % sizeof(int)) == 0);
    if (fs == 0) return;
    size_t ns = fs / sizeof(int);
    int* cns_seq_id_array = new int[ns];
    hbn_dfopen(in, path, "rb");
    hbn_fread(cns_seq_id_array, sizeof(int), ns, in);
    hbn_fclose(in);
    in = NULL;

    int lid = *local_id;
    for (int i = 0; i < ns; ++i) {
        int id = cns_seq_id_array[i];
        hbn_assert(id_map_array[id].first == id);
        hbn_assert(id_map_array[id].second == -1);
        id_map_array[id].second = lid++;
    }
    delete[] cns_seq_id_array;
    *local_id = lid;
}

static pair<int, int>*
s_make_id_map(const char* old_db_dir, const char* old_cns_sr_dir)
{
    int old_num_reads = seqdb_load_num_reads(old_db_dir, INIT_QUERY_DB_TITLE);
    pair<int, int>* id_map_array = new pair<int, int>[old_num_reads];
    for (int i = 0; i < old_num_reads; ++i) {
        id_map_array[i].first = i;
        id_map_array[i].second = -1;
    }
    int old_np = load_partition_count(old_cns_sr_dir, NULL);
    int local_id = 0;
    for (int pid = 0; pid < old_np; ++pid) {
        s_add_one_part_seq_ids(old_cns_sr_dir, pid, &local_id, id_map_array);
    }
    return id_map_array;
}

static void
s_map_one_sr_part(const char* old_sr_dir,
    const char* new_sr_dir,
    const int pid,
    pair<int, int>* id_map_array)
{
    char job_name[256];
    char pid_str[64];
    u64_to_fixed_width_string_r(pid, pid_str, HBN_DIGIT_WIDTH);
    sprintf(job_name, "make_new_sr_%s", pid_str);
    hbn_timing_begin(job_name);
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(old_sr_dir, DEFAULT_PART_PREFIX, pid, path);
    size_t hit_count = 0;
    HbnConsensusInitHit* hit_array = (HbnConsensusInitHit*)load_part_records(path, sizeof(HbnConsensusInitHit), &hit_count);
    make_partition_name(new_sr_dir, DEFAULT_PART_PREFIX, pid, path);
    hbn_dfopen(out, path, "wb");
    for (size_t i = 0; i < hit_count; ++i) {
        HbnConsensusInitHit hit = hit_array[i];
        int qid = hit.qid;
        int sid = hit.sid;
        if (id_map_array[qid].second == -1 || id_map_array[sid].second == -1) continue;
        qid = id_map_array[qid].second;
        sid = id_map_array[sid].second;
        hit.qid = qid;
        hit.sid = sid;
        hbn_fwrite(&hit, sizeof(HbnConsensusInitHit), 1, out);
    }
    hbn_fclose(out);
    free(hit_array);
    hbn_timing_end(job_name);
}

static void
print_usage(const char* pn)
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s old_db_dir old_sr_dir new_sr_di\n", pn);
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }
    const char* old_db_dir = argv[1];
    const char* old_sr_dir = argv[2];
    const char* new_sr_dir = argv[3];
    create_directory(new_sr_dir);
    pair<int, int>* id_map_array = s_make_id_map(old_db_dir, old_sr_dir);
    int np = load_partition_count(old_sr_dir, NULL);
    dump_partition_count(new_sr_dir, NULL, np);
    for (int i = 0; i < np; ++i) {
        s_map_one_sr_part(old_sr_dir, new_sr_dir, i, id_map_array);
    }
    delete[] id_map_array;
    return 0;
}