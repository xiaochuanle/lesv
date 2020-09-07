#include "../../corelib/m4_record.h"
#include "../../corelib/name2id_map.h"
#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <fstream>
#include <string>
#include <cassert>

using namespace std;

void
read_one_m4(NameToIdMap* query_name2id, NameToIdMap* subject_name2id, const string& line, M4Record* m4)
{
    NStr::CTempString kDelim("\t");
    vector<string> components;
    NStr::Split(line, kDelim, components);
    hbn_assert(components.size() == 12);
    m4->qid = Name2IdMap_name2id(query_name2id, components[0].c_str());
    m4->sid = Name2IdMap_name2id(subject_name2id, components[1].c_str());
    m4->ident_perc = NStr::StringToDouble(components[2]);
    m4->score = NStr::StringToInt(components[3]);
    m4->qdir = NStr::StringToInt(components[4]);
    m4->qoff = NStr::StringToInt(components[5]);
    m4->qend = NStr::StringToInt(components[6]);
    m4->qsize = NStr::StringToInt(components[7]);
    m4->sdir = NStr::StringToInt(components[8]);
    m4->soff = NStr::StringToInt(components[9]);
    m4->send = NStr::StringToInt(components[10]);
    m4->ssize = NStr::StringToInt(components[11]);
}

void
m4x_one_volume(NameToIdMap* query_name2id, NameToIdMap* subject_name2id, const char* db_dir, const int qi, const int sj)
{
    char path[HBN_MAX_PATH_LEN];
    char qi_buf[64];
    char sj_buf[64];
    u64_to_fixed_width_string_r(qi, qi_buf, HBN_DIGIT_WIDTH);
    u64_to_fixed_width_string_r(sj, sj_buf, HBN_DIGIT_WIDTH);
    sprintf(path, "%s/backup_results/stagebackup_results_Q%s_D%s", db_dir, qi_buf, sj_buf);
    char job_name[256];
    sprintf(job_name, "Q%s_D%s", qi_buf, sj_buf);
    hbn_timing_begin(job_name);
    HBN_LOG("%s", path);
    ifstream in(path);
    assert(in.is_open());
    strcat(path, ".m4x");
    hbn_dfopen(out, path, "wb");
    string line;
    M4Record m4;
    int nm4 = 0;
    while (getline(in, line)) {
        read_one_m4(query_name2id, subject_name2id, line, &m4);
        hbn_fwrite(&m4, sizeof(M4Record), 1, out);
        ++nm4;
    }
    hbn_fclose(out);
    HBN_LOG("number of m4: %d", nm4);
    hbn_timing_end(job_name);
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s db_dir\n", argv[0]);
        return 1;
    }
    const char* db_dir = argv[1];
    CSeqDBInfo query_dbinfo = seqdb_load_volume_info(db_dir, INIT_QUERY_DB_TITLE, 0);
    CSeqDBInfo subject_dbinfo = seqdb_load_volume_info(db_dir, INIT_SUBJECT_DB_TITLE, 0);
    char* query_names = load_seq_headers(db_dir, INIT_QUERY_DB_TITLE, query_dbinfo.hdr_offset_from, query_dbinfo.hdr_offset_to);
    char* subject_names = load_seq_headers(db_dir, INIT_SUBJECT_DB_TITLE, subject_dbinfo.hdr_offset_from, subject_dbinfo.hdr_offset_to);
    NameToIdMap* query_name2id = NameToIdMapNew();
    NameToIdMapSet(query_names, query_dbinfo.num_seqs, query_name2id);
    NameToIdMap* subject_name2id = NameToIdMapNew();
    NameToIdMapSet(subject_names, subject_dbinfo.num_seqs, subject_name2id);

    int num_query_vols = seqdb_load_num_volumes(db_dir, INIT_QUERY_DB_TITLE);
    int num_subject_vols = seqdb_load_num_volumes(db_dir, INIT_SUBJECT_DB_TITLE);
    for (int i = 0; i < num_query_vols; ++i) {
        for (int j = 0; j < num_subject_vols; ++j) {
            m4x_one_volume(query_name2id, subject_name2id, db_dir, i, j);
        }
    }

    NameToIdMapFree(subject_name2id);
    NameToIdMapFree(query_name2id);
    free(subject_names);
    free(query_names);
}
