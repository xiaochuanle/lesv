#include "../../corelib/build_db.h"
#include "../../corelib/cmd_arg.h"
#include "../../corelib/seqdb.h"
#include "../../corelib/hbn_package_version.h"

#include <sys/stat.h>

static void
dump_usage(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s db_dir db_title\n", pn);
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        dump_usage(argv[0]);
        return EXIT_FAILURE;
    }
    const char* db_dir = argv[1];
    const char* db_title = argv[2];
    const int num_vols = seqdb_load_num_volumes(db_dir, db_title);
    kv_dinit(vec_u8, seq_v);
    for (int i = 0; i < num_vols; ++i) {
        HBN_LOG("dumpping volume %d", i);
        CSeqDB* vol = seqdb_load(db_dir, db_title, i);
        for (int j = 0; j < vol->dbinfo.num_seqs; ++j) {
            const char* name = seqdb_seq_name(vol, j);
            seqdb_extract_sequence(vol, j, FWD, &seq_v);
            for (size_t p = 0; p < kv_size(seq_v); ++p) {
                int c = kv_A(seq_v, p);
                c = DECODE_RESIDUE(c);
                kv_A(seq_v, p) = c;
            }
            fprintf(stdout, ">%s\n", name);
            hbn_fwrite(kv_data(seq_v), 1, kv_size(seq_v), stdout);
            fprintf(stdout, "\n");
        }
        CSeqDBFree(vol);
    }
    kv_destroy(seq_v);
}