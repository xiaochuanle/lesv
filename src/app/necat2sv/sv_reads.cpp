#include "sv_reads.h"

#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <fstream>

using namespace std;

extern "C"
void
sread_sv_read(const char* in, SvRead* rp)
{
    ncbi::CTempString ins(in);
    ncbi::CTempString delim("\t");
    vector<ncbi::CTempString> components;
    NStr::Split(ins, delim, components);
    int np = components.size();
    hbn_assert(np == 9 || np == 10 || np == 15 || np == 16, "n = %d, %s", np, in);
    rp->query_id = NStr::StringToInt(components[0]);
    rp->qdir = NStr::StringToInt(components[1]);
    rp->qoff = NStr::StringToInt(components[2]);
    rp->qend = NStr::StringToInt(components[3]);
    rp->qsize = NStr::StringToInt(components[4]);
    rp->subject_id = NStr::StringToInt(components[5]);
    rp->soff = NStr::StringToInt(components[6]);
    rp->send = NStr::StringToInt(components[7]);
    rp->dist = NStr::StringToInt(components[8]);

    rp->dual_qdir = -1;
    if (np > 11) {
        rp->dual_qdir = NStr::StringToInt(components[9]);
        rp->dual_qoff = NStr::StringToInt(components[10]);
        rp->dual_qend = NStr::StringToInt(components[11]);
        rp->dual_soff = NStr::StringToInt(components[12]);
        rp->dual_send = NStr::StringToInt(components[13]);
        rp->dual_dist = NStr::StringToInt(components[14]);
    }
}

extern "C"
SvRead*
load_sv_read_array(const char* path, int* sv_read_count)
{
    ifstream in(path);
    if (!in) HBN_ERR("Fail to open file %s for reading", path);
    string line;
    kv_dinit(vec_sv_read, sv_read_list);
    SvRead sv_read;
    while (getline(in, line)) {
        sread_sv_read(line.c_str(), &sv_read);
        kv_push(SvRead, sv_read_list, sv_read);
    }
    in.close();
    *sv_read_count = kv_size(sv_read_list);
    return kv_data(sv_read_list);
}