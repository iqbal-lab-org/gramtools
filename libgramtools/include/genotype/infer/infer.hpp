/**
 * @file
 * Interface to running personalised reference inference and
 * genotyping from coverage-annotated PRG.
 */

#ifndef GRAMTOOLS_INFER_HPP
#define GRAMTOOLS_INFER_HPP

#include "genotype/infer/genotyping_models.hpp"
#include "common/read_stats.hpp"

using namespace gram::genotype::infer;
using namespace gram::genotype::infer::probabilities;

class LevelGenotyper{
   coverage_Graph const* cov_graph;
   SitesGroupedAlleleCounts const* gped_covs;
   likelihood_related_stats* l_stats;

   PoissonLogPmf poisson_prob;
   gt_sites genotyped_records;

public:
    LevelGenotyper() : cov_graph(nullptr), gped_covs(nullptr) {}
    LevelGenotyper(coverage_Graph const& cov_graph, SitesGroupedAlleleCounts const& gped_covs,
                   ReadStats const& read_stats);
    CovCount find_minimum_non_error_cov(double mean_pb_error, poisson_pmf_ptr poisson_prob);
};

#endif //GRAMTOOLS_INFER_HPP
