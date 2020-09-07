#ifndef __LOF_HPP
#define __LOF_HPP

#include "lof_point.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include <cassert>
#include <cfloat>
#include <iostream>

namespace hbn_lof {

using namespace std;

#define LOG_BK_INFO "[" << __FILE__ << ", " << __FUNCTION__ << ", " << __LINE__ << "] "

static bool
s_double_equal(double x, double y)
{
    return fabs(x - y) < 1.0e-6;
}

template <int DIM>
double distance_euclidean(const LofPoint<DIM>& x, const LofPoint<DIM>& y)
{
    double sum = 0.0;
    for (int i = 0; i < DIM; ++i) {
        double xi = x[i];
        double yi = y[i];
        sum += (xi - yi) * (xi - yi);
    }
    return sqrt(sum/DIM);
}

template <int DIM>
void compute_instance_attribute_bounds(const vector<LofPoint<DIM>>& instances,
        array<double, DIM>& max_attribute_values,
        array<double, DIM>& min_attribute_values)
{
    fill(max_attribute_values.begin(), max_attribute_values.end(), numeric_limits<double>::min());
    fill(min_attribute_values.begin(), min_attribute_values.end(), numeric_limits<double>::max());
    for (const auto& inst : instances) {
        for (int i = 0; i < DIM; ++i) {
            double xi = inst[i];
            max_attribute_values[i] = max(max_attribute_values[i], xi);
            min_attribute_values[i] = min(min_attribute_values[i], xi);
        }
    }
}

template <int DIM>
void normalise_instance(LofPoint<DIM>& instance,
        array<double, DIM>& max_attribute_values,
        array<double, DIM>& min_attribute_values)
{
    for (int i = 0; i < DIM; ++i) {
        double max_xi = max_attribute_values[i];
        double min_xi = min_attribute_values[i];
        double n_xi = .0;
        if (s_double_equal(max_xi, min_xi)) {
            n_xi = .0;
        } else {
            double xi = instance[i];
            n_xi = (xi - min_xi) / (max_xi - min_xi);
        }
        instance[i] = n_xi;
    }
}

template <int DIM>
void normalise_instances(vector<LofPoint<DIM>>& instances,
        array<double, DIM>& max_attribute_values,
        array<double, DIM>& min_attribute_values)
{
    for (auto& inst : instances) normalise_instance<DIM>(inst, max_attribute_values, min_attribute_values);
}

template <int DIM>
double k_distance(const int k, 
            LofPoint<DIM>& instance,
            vector<LofPoint<DIM>>& instances,
            vector<LofPoint<DIM>>& nbhd_instances)
{
    vector<pair<int, double>> idx_and_dist_list;
    int i = 0;
    //cerr << "target instance: " << instance << endl;
    for (auto& instance2 : instances) {
        double distance_value = distance_euclidean(instance, instance2);
        //cerr << instance2 << ", euclid_dist = " << distance_value << endl;
        idx_and_dist_list.push_back(pair<int, double>(i, distance_value));
        ++i;
    }
    sort(idx_and_dist_list.begin(),
        idx_and_dist_list.end(),
        [](const pair<int, double>& x, const pair<int, double>& y)->bool { return x.second < y.second; });

    const int n = idx_and_dist_list.size();
    double k_distance_value = 0.0;
    int k_i = 0;
    i = 0;
    nbhd_instances.clear();
    while (i < n) {
        double di = idx_and_dist_list[i].second;
        int j = i + 1;
        while (j < n) {
            double dj = idx_and_dist_list[j].second;
            if (!s_double_equal(di, dj)) break;
            ++j;
        }
        for (int p = i; p < j; ++p) {
            int x = idx_and_dist_list[p].first;
            nbhd_instances.push_back(instances[x]);
        }
        ++k_i;
        k_distance_value = di;
        if (k_i == k) break;
        i = j;
    }
    return k_distance_value;
}

template <int DIM>
double reachability_distance(const int k,
            LofPoint<DIM>& instance1,
            LofPoint<DIM>& instance2,
            vector<LofPoint<DIM>>& instances)
{
    vector<LofPoint<DIM>> nbhd_instances;
    double k_distance_value = k_distance(k, instance2, instances, nbhd_instances);
    double euclid_distance = distance_euclidean(instance1, instance2);
    return max(k_distance_value, euclid_distance);
}

template <int DIM>
double local_reachability_density(const int min_ptrs,
            LofPoint<DIM>& instance,
            vector<LofPoint<DIM>>& instances)
{
    vector<LofPoint<DIM>> nbhd_instances;
    double k_distance_value = k_distance(min_ptrs, instance, instances, nbhd_instances);
    //cerr << instance << '\t' << k_distance_value << endl;
    double sum_reachability_distance = 0.0;
    for (auto& neighbour : nbhd_instances) {
        sum_reachability_distance += reachability_distance(min_ptrs,
            instance, 
            neighbour,
            instances);
    }
    //cerr << "number of instances: " << instances.size() << endl;
    //cerr << instance << '\t';
    //cerr << "k_dist = " << k_distance_value << ", sum_reach_dist = " << sum_reachability_distance << endl;
    if (s_double_equal(sum_reachability_distance, 0.0)) {
        cerr << LOG_BK_INFO 
             << "Instance " << instance
             << " (could be normalized) is identical to all the neighbors."
             << " Setting local reachability density to inf"
             << endl;
        return DBL_MIN;;
    }  else {
        return static_cast<double>(nbhd_instances.size()) / sum_reachability_distance;
    }
}

template <int DIM>
double local_outlier_factor(const int min_pts,
            LofPoint<DIM>& instance,
            vector<LofPoint<DIM>>& instances)
{
    vector<LofPoint<DIM>> nbhd_instances;
    vector<LofPoint<DIM>> instances_without_nbhd;
    double k_distance_value = k_distance(min_pts, instance, instances, nbhd_instances);
    //cerr << "number of neighbour: " << nbhd_instances.size() << endl;
    double instance_lrd = local_reachability_density(min_pts, instance, instances);
    //cerr << "k_dist = " << k_distance_value << ", lrd = " << instance_lrd << endl;
    double sum_lrd_ratios = 0.0;
    //cerr << "total instances: " << instances.size() << endl;
    for (size_t i = 0; i < nbhd_instances.size(); ++i) {
        instances_without_nbhd.clear();
        LofPoint<DIM>& neighbour = nbhd_instances[i];
        for (auto& inst : instances) {
            if (neighbour != inst) {
                instances_without_nbhd.push_back(inst);
            }
        }
        double neighbour_lrd = local_reachability_density(min_pts, neighbour, instances_without_nbhd);
        sum_lrd_ratios += neighbour_lrd / instance_lrd;
        //cerr << "nbhd " << i << '\t' << neighbour << endl;
        //cerr << "lrd = " << neighbour_lrd << endl;
        //cerr << "number of instances: " << instances_without_nbhd.size() << endl;        
    }
    //cerr << "k_dist = " << k_distance_value << ", sum_lrd_ratios = " << sum_lrd_ratios << endl;
    return sum_lrd_ratios / static_cast<double>(nbhd_instances.size());
}

template <int DIM>
double lof_local_outlier_factor(const int min_pts,
            LofPoint<DIM>& instance,
            vector<LofPoint<DIM>>& instances,
            const bool do_normalise)
{
    if (do_normalise) {
        array<double, DIM> max_attribute_values;
        array<double, DIM> min_attribute_values;
        compute_instance_attribute_bounds<DIM>(instances, max_attribute_values, min_attribute_values);
        normalise_instances<DIM>(instances, max_attribute_values, min_attribute_values);
        normalise_instance<DIM>(instance, max_attribute_values, min_attribute_values);
    }
    return local_outlier_factor<DIM>(min_pts, instance, instances);
}

template <int DIM>
struct LofOutlier
{
    double lof;
    int index;
    LofPoint<DIM> instance;

    LofOutlier(double lof_, int index_, LofPoint<DIM>& instance_)
        : lof(lof_), index(index_), instance(instance_) {}
};

template <int DIM>
vector<LofOutlier<DIM>>
GetOutliers(int min_pts, vector<LofPoint<DIM>>& instances)
{
    if (0)
    {
        array<double, DIM> max_coords;
        array<double, DIM> min_coords;
        compute_instance_attribute_bounds<DIM>(instances, max_coords, min_coords);
        normalise_instances<DIM>(instances, max_coords, min_coords);
    }

    vector<LofOutlier<DIM>> outliers;
    vector<LofPoint<DIM>> backup_instances;
    for (size_t i = 0; i < instances.size(); ++i) {
        backup_instances = instances;
        LofPoint<DIM> instance = backup_instances[i];
        //backup_instances[i] = backup_instances.back();
        //backup_instances.pop_back();
        backup_instances.erase(backup_instances.begin() + i);
        //cerr << i << '\t' << instance;
        double value = lof_local_outlier_factor(min_pts, instance, backup_instances, false);
        //cerr << '\t' << value << endl;
        if (value > 1.0) {
            outliers.push_back(LofOutlier<DIM>(value, i, instances[i]));
        }
        //break;
    }
    sort(outliers.begin(),
        outliers.end(),
        [](const LofOutlier<DIM>& x, const LofOutlier<DIM>& y)->bool { return x.lof > y.lof; });
    return outliers;
}

} // hbn_lof

#endif // __LOF_HPP
