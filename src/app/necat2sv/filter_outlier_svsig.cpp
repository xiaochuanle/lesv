#include "filter_outlier_svsig.h"

#include "../../algo/local_outlier_factor/lof.hpp"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include <iostream>
#include <set>
#include <map>

using namespace hbn_lof;
using namespace std;

extern "C"
void
find_outlier_svsig(SvSignature* svsig_array, int svsig_count)
{
#if 1
    vector<LofPoint<2>> pos_and_len_list;
    for (int i = 0; i < svsig_count; ++i) {
        int len = svsig_array[i].qto - svsig_array[i].qfrom;
        int soff = svsig_array[i].sfrom - svsig_array[0].sfrom;
        LofPoint<2> pnt(soff, len);
        pos_and_len_list.push_back(pnt);
    }   
    vector<LofOutlier<2>> outliers = GetOutliers<2>(5, pos_and_len_list);
    HBN_LOG("total svsigs: %d, number of outliers: %zu", svsig_count, outliers.size());
    for (auto& outlier : outliers) {
        int i = outlier.index;
        dump_svsig(fprintf, stderr, svsig_array[i], NULL);
        int len = svsig_array[i].qto - svsig_array[i].qfrom;
        fprintf(stderr, "len = %d\n", len);
        svsig_array[i].qid = -1;
    }
    HBN_LOG("normal svsigs:");
    for (int i = 0; i < svsig_count; ++i) {
        if (svsig_array[i].qid == -1) continue;
        dump_svsig(fprintf, stderr, svsig_array[i], NULL);
        int len = svsig_array[i].qto - svsig_array[i].qfrom;
        fprintf(stderr, "len = %d\n", len);
    }
#endif 

#if 0
    vector<LofPoint<1>> len_list;
    for (int i = 0; i < svsig_count; ++i) {
        int len = svsig_array[i].qto - svsig_array[i].qfrom;
        LofPoint<1> pnt(len);
        len_list.push_back(pnt);
    }
    vector<LofOutlier<1>> outliers = GetOutliers<1>(5, len_list);
    HBN_LOG("number of outliers: ", outliers.size());
    for (auto& outlier : outliers) {
        int i = outlier.index;
        dump_svsig(fprintf, stderr, svsig_array[i], NULL);
        int len = svsig_array[i].qto - svsig_array[i].qfrom;
        fprintf(stderr, "len = %d\n", len);
        svsig_array[i].qid = -1;
    }
    HBN_LOG("normal svsigs:");
    for (int i = 0; i < svsig_count; ++i) {
        if (svsig_array[i].qid == -1) continue;
        dump_svsig(fprintf, stderr, svsig_array[i], NULL);
        int len = svsig_array[i].qto - svsig_array[i].qfrom;
        fprintf(stderr, "len = %d\n", len);
    }    
#endif
}

int calc_avg_cov(int* cov_stats, int from, int to)
{
    int sum = 0;
    for (int i = from; i < to; ++i) sum += cov_stats[i];
    return sum / (to - from);
}

extern "C"
void trim_svsig_by_cov_stats(SvSignature* svga, int svgc)
{
    int size = svga[0].ssize;
    int* cov_stats = new int[size];
    fill(cov_stats, cov_stats + size, 0);
    set<int> added_ids;
    size_t total_res = 0;
    for (int i = 0; i < svgc; ++i) {
        if (added_ids.find(svga[i].qid) != added_ids.end()) continue;
        for (int p = svga[i].fsfrom; p < svga[i].fsto; ++p) {
            ++cov_stats[p];
        }
        total_res += (svga[i].fqto - svga[i].fqfrom);
        added_ids.insert(svga[i].qid);
    }
    string total_res_str = NStr::UInt8ToString_DataSize(total_res);
    HBN_LOG("%zu distinct reads (%s)", added_ids.size(), total_res_str.c_str());
    int i = 0;
    const int kMaxCov = 80;
    while (i < size) {
        while (i < size && cov_stats[i] < kMaxCov) ++i;
        if (i >= size) break;
        int j = i + 1;
        while (j < size && cov_stats[j] >= kMaxCov) ++j;
        int avg_cov = calc_avg_cov(cov_stats, i, j);
        HBN_LOG("%d --- %d, avg_cov = %d", i, j, avg_cov);
        i = j;
    }
}