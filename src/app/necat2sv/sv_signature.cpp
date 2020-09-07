#include "sv_signature.h"

#include <boost/algorithm/sorting/spread_sort.hpp>
#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <fstream>

using namespace boost;
using namespace std;

/////////////////

struct SvSignature_SfromLessThan {
    inline bool operator()(const SvSignature &x, const SvSignature &y) const { return x.sfrom < y.sfrom; }
};

struct SvSignature_SfromRightShift {
        inline int operator()(const SvSignature &x, const unsigned offset) { return x.sfrom >> offset; }
};

extern "C"
void
sort_svsig_sfrom_lt(size_t n, SvSignature* a)
{
    integer_sort(a, a + n, SvSignature_SfromRightShift(), SvSignature_SfromLessThan());
}

extern "C"
void
sread_sv_signature(const char* in, SvSignature* svsig)
{
    ncbi::CTempString ins(in);
    ncbi::CTempString delim("\t");
    vector<ncbi::CTempString> components;
    NStr::Split(ins, delim, components);
    int np = components.size();
    hbn_assert(np == 13 || np == 14);
    svsig->qid = NStr::StringToInt(components[0]);
    svsig->qdir = NStr::StringToInt(components[1]);
    svsig->qfrom = NStr::StringToInt(components[2]);
    svsig->qto = NStr::StringToInt(components[3]);
    svsig->fqfrom = NStr::StringToInt(components[4]);
    svsig->fqto = NStr::StringToInt(components[5]);
    svsig->qsize = NStr::StringToInt(components[6]);
    svsig->sfrom = NStr::StringToInt(components[7]);
    svsig->sto = NStr::StringToInt(components[8]);
    svsig->fsfrom = NStr::StringToInt(components[9]);
    svsig->fsto = NStr::StringToInt(components[10]);
    svsig->ssize = NStr::StringToInt(components[11]);
    svsig->type = static_cast<EGapAlignOpType>(NStr::StringToInt(components[12]));
    hbn_assert(svsig->type >= eGapAlignDel && svsig->type < eGapAlignInvalid);
}

extern "C"
SvSignature*
load_sv_signature_array(const char* path, int* sv_signature_count, const EGapAlignOpType type)
{
    ifstream in(path);
    hbn_assert(in);
    string line;
    kv_dinit(vec_sv_sig, sv_sig_list);
    SvSignature svsig;
    while (getline(in, line)) {
        sread_sv_signature(line.c_str(), &svsig);
        if (type == eGapAlignInvalid || svsig.type == type) {
            kv_push(SvSignature, sv_sig_list, svsig);
        }
    }
    in.close();
    *sv_signature_count = kv_size(sv_sig_list);
    return kv_data(sv_sig_list);    
}
