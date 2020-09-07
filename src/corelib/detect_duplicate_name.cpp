#include "detect_duplicate_name.h"

#include <map>

using namespace std;

struct CmpCStr
{
    inline bool operator()(const char* a, const char* b) const {
        return strcmp(a, b) < 0;
    }
};

extern "C"
void detect_duplicate_sequence_names(const char* db_dir, const char* db_title, const char* fasta_input_path)
{
    CSeqDBInfo dbinfo = seqdb_load_volume_info(db_dir, db_title, 0);
    hbn_assert(dbinfo.seq_start_id == 0);
    if (dbinfo.num_seqs == 0) return;

    CSeqInfo* seqinfo_array = load_seq_infos(db_dir, db_title, dbinfo.seq_start_id, dbinfo.num_seqs);
    char* seqhdr_array = load_seq_headers(db_dir, db_title, dbinfo.hdr_offset_from, dbinfo.hdr_offset_to);

    map<const char*, int, CmpCStr> added_names;
    for (int i = dbinfo.seq_start_id; i < dbinfo.seq_start_id + dbinfo.num_seqs; ++i) {
        const char* name = seqhdr_array + seqinfo_array[i].hdr_offset;
        auto name_pos = added_names.find(name);
        if (name_pos != added_names.end()) {
            if ((*name_pos).second == 0) {
                fprintf(stderr, "[HbnFastaReader] Warning: Sequence name '%s' is duplicate in DB '%s' (If a name occurs more than two times, it will be warned only once.)\n",
                    name, fasta_input_path);
                (*name_pos).second = 1;
            }
        } else {
            added_names.insert(pair<const char*, int>(name, 0));
        }
    }

    free(seqinfo_array);
    free(seqhdr_array);
}