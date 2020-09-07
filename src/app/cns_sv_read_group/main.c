#include "cmdline_args.h"
#include "cns_aux.h"
#include "cns_one_group.h"
#include "../necat2sv/sv_read_group_file_name.h"

#include <errno.h>
#include <pthread.h>
#include <sys/stat.h>

int main(int argc, char* argv[])
{
    HbnProgramOptions opts;
    ParseHbnProgramCmdLineArguments(argc, argv, &opts);
    CSeqDBInfo raw_read_dbinfo = seqdb_load_volume_info(opts.db_dir, INIT_QUERY_DB_TITLE, 0);
    CSeqInfo* raw_read_info_array = load_seq_infos(opts.db_dir, INIT_QUERY_DB_TITLE, 0, raw_read_dbinfo.num_seqs);
    char* raw_read_name_array = load_seq_headers(opts.db_dir, INIT_QUERY_DB_TITLE, 0, raw_read_dbinfo.hdr_offset_to);
    char path[HBN_MAX_PATH_LEN];
    make_packed_seq_path(opts.db_dir, INIT_QUERY_DB_TITLE, path);
    hbn_dfopen(pac_file, path, "rb");
    CSeqDB* db = seqdb_load(opts.db_dir, INIT_SUBJECT_DB_TITLE, 0);
    for (int i = 0; i < db->dbinfo.num_seqs; ++i) {
        //if (i != 20) continue;
        if (subject_sv_read_group_is_corrected(opts.wrk_dir, i)) {
            HBN_LOG("sv read groups for subject %d are corrected.", i);
            continue;
        }
        cns_suject_svr_groups_st(argc, argv, 
            &opts, db, raw_read_info_array, 
            raw_read_name_array, pac_file, i);
        subject_sv_read_group_make_corrected(opts.wrk_dir, i);
    }
    hbn_fclose(pac_file);
    sfree(raw_read_name_array);
    sfree(raw_read_info_array);
    CSeqDBFree(db);
    return 0;
}
