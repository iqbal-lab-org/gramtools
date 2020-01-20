#include <cmath>
#include "genotype/infer/infer.hpp"

using namespace gram::genotype::infer::probabilities;

CovCount LevelGenotyper::find_minimum_non_error_cov(double mean_pb_error, poisson_pmf_ptr poisson_prob) {
    CovCount min_count{1};
    while ( (*poisson_prob)(params{static_cast<double>(min_count)}) <= log(pow(mean_pb_error, min_count)) )
        ++min_count;
    return min_count;
}

LevelGenotyper::LevelGenotyper(coverage_Graph const& cov_graph, SitesGroupedAlleleCounts const& gped_covs,
                               ReadStats const& read_stats) :
        cov_graph(&cov_graph), gped_covs(&gped_covs){
   auto mean_cov_depth = read_stats.get_mean_cov_depth();
   auto mean_pb_error = read_stats.get_mean_pb_error();
   poisson_prob = PoissonLogPmf{params{mean_cov_depth}};

   likelihood_related_stats new_l_stats{
       mean_cov_depth, mean_pb_error,
       log(1 - exp(-mean_cov_depth)),
       log(1 - exp(-mean_cov_depth / 2)),
       find_minimum_non_error_cov(mean_pb_error, &poisson_prob)
   };
}
