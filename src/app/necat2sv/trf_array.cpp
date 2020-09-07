#include "trf_array.h"

#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include <string>
#include <fstream>

using namespace std;

TrfArray*
TrfArrayBuild(const char* trf_path, const char* db_dir, const CSeqInfo* seqinfo_array)
{
    TrfArray* array = (TrfArray*)calloc(1, sizeof(TrfArray));
    array->dbinfo = seqdb_load_volume_info(db_dir, INIT_SUBJECT_DB_TITLE, 0);
    array->seqinfo_array = seqinfo_array;
    const size_t max_res_offset = seqinfo_array[array->dbinfo.num_seqs-1].seq_offset
                                  +
                                  seqinfo_array[array->dbinfo.num_seqs-1].seq_size;
    array->trf_array = (char*)calloc(sizeof(char), max_res_offset);

    char* subject_names = load_seq_headers(db_dir, 
                            INIT_SUBJECT_DB_TITLE, 
                            array->dbinfo.hdr_offset_from, 
                            array->dbinfo.hdr_offset_to);
    NameToIdMap* name2id_map = NameToIdMapNew();
    NameToIdMapSet(subject_names, array->dbinfo.num_seqs, name2id_map);

    ifstream in(trf_path);
    hbn_assert(in);
    string line;
    NStr::CTempString kDelim("\t");
    vector<string> components;
    size_t trf_res = 0;
    while (getline(in, line)) {
        components.clear();
        NStr::Split(line, kDelim, components);
        hbn_assert(components.size() == 3);
        if (!Name2IdMap_NameExists(name2id_map, components[0].c_str())) continue;
        int sid = Name2IdMap_name2id(name2id_map, components[0].c_str());
        hbn_assert(sid >= 0);
        hbn_assert(sid < array->dbinfo.num_seqs);
        size_t soff = seqinfo_array[sid].seq_offset;
        size_t from = NStr::StringToUInt8(components[1]);
        size_t to = NStr::StringToUInt8(components[2]);
        hbn_assert(from < to);
        hbn_assert(to <= seqinfo_array[sid].seq_size);
        trf_res += (to - from);
        from += soff;
        to += soff;
        for (size_t i = from; i < to; ++i) array->trf_array[i] = 1;
    }
    NameToIdMapFree(name2id_map);
    free(subject_names);
    
    size_t total_res = 0;
    for (int i = 0; i < array->dbinfo.num_seqs; ++i) total_res += seqinfo_array[i].seq_size;
    string total_res_str = NStr::UInt8ToString_DataSize(total_res);
    string trf_res_str = NStr::UInt8ToString_DataSize(trf_res);
    double trf_perc = 100.0 * trf_res / total_res;
    string trf_perc_str = NStr::DoubleToString(trf_perc);
    HBN_LOG("Total residues: %s, trf residues: %s (%s%%)", 
        total_res_str.c_str(), trf_res_str.c_str(), trf_perc_str.c_str());

    return array;
}

TrfArray*
TrfArrayFree(TrfArray* array)
{
    free(array->trf_array);
    free(array);
    return NULL;
}

BOOL
fall_in_trf_subseq(TrfArray* array, int sid, size_t from, size_t to)
{
    hbn_assert(sid < array->dbinfo.num_seqs);
    size_t soff = array->seqinfo_array[sid].seq_offset;
    from += soff;
    to += soff;
    size_t trf_cnt = 0;
    for (size_t i = from; i < to; ++i) trf_cnt += array->trf_array[i];
    return (to - from - trf_cnt <= 2000);
}