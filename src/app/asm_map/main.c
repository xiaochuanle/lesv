#include "cmdline_args.h"
#include "map_one_part.h"
#include "../../corelib/raw_reads_reader.h"
#include "../../corelib/partition_mt.h"

#include <errno.h>
#include <sys/stat.h>

int main(int argc, char* argv[])
{
    HbnProgramOptions opts;
    ParseHbnProgramCmdLineArguments(argc, argv, &opts);
    int num_parts = load_partition_count(opts.sr_dir, NULL);
    RawReadReader* raw_reads = RawReadReaderNew(opts.db_dir, INIT_QUERY_DB_TITLE, opts.use_batch_mode);
    char job_name[256];
    char pid_str[64];
    for (int i = opts.node_id; i < num_parts; i += opts.num_nodes) {
        u64_to_fixed_width_string_r(i, pid_str, HBN_DIGIT_WIDTH);
        sprintf(job_name, "mapping part %s", pid_str);
        hbn_timing_begin(job_name);
        asm_map_one_part(opts.sr_dir, &opts, raw_reads, i, opts.batch_size, opts.use_batch_mode);
        hbn_timing_end(job_name);
    }
    raw_reads = RawReadReaderFree(raw_reads);
    return 0;
}