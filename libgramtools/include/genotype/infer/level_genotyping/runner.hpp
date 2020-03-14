/**
 * @file
 * Interface to running personalised reference inference and
 * genotyping from coverage-annotated PRG.
 *
 * Although a PRG moves beyond a single reference genome, each site's first allele will consistently
 * be the first allele (haplogroup) of the bubble in the graph, so that we can use that as REF if needed.
 */

#ifndef LVLGT_RUNNER
#define LVLGT_RUNNER

#include "model.hpp"
#include "genotype/read_stats.hpp"

using namespace gram::genotype::infer;
using namespace gram::genotype::infer::probabilities;

namespace gram::genotype::infer {


class LevelGenotyper : public Genotyper {
    likelihood_related_stats l_stats;
    Ploidy ploidy;

public:
    LevelGenotyper() = default;
    LevelGenotyper(child_map const& ch, gt_sites const& sites) : Genotyper(sites, ch) {}

    LevelGenotyper(coverage_Graph const &cov_graph, SitesGroupedAlleleCounts const &gped_covs,
                   ReadStats const &read_stats, Ploidy const ploidy);

    header_vec get_model_specific_headers() override;

    static CovCount find_minimum_non_error_cov(double mean_pb_error, poisson_pmf_ptr poisson_prob);
    static likelihood_related_stats make_l_stats(double mean_cov_depth, double mean_pb_error);

};
}

#endif //LVLGT_RUNNER
