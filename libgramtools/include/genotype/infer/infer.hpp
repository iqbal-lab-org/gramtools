/**
 * @file
 * Interface to running personalised reference inference and
 * genotyping from coverage-annotated PRG.
 *
 * Although a PRG moves beyond a single reference genome, each site's first allele will consistently
 * be the first allele (haplogroup) of the bubble in the graph, so that we can use that as REF if needed.
 */

#ifndef GRAMTOOLS_INFER_HPP
#define GRAMTOOLS_INFER_HPP

#include "genotype/infer/genotyping_models.hpp"
#include "genotype/read_stats.hpp"

using namespace gram::genotype::infer;
using namespace gram::genotype::infer::probabilities;

namespace gram::genotype::infer {


class LevelGenotyper {
    coverage_Graph const *cov_graph;
    SitesGroupedAlleleCounts const *gped_covs;
    likelihood_related_stats l_stats;

    Ploidy ploidy;
    gt_sites genotyped_records;
    child_map child_m;

public:
    LevelGenotyper() : cov_graph(nullptr), gped_covs(nullptr) {}
    LevelGenotyper(child_map const& ch, gt_sites const& sites) :
        child_m(ch), genotyped_records(sites), cov_graph(nullptr), gped_covs(nullptr) {}

    LevelGenotyper(coverage_Graph const &cov_graph, SitesGroupedAlleleCounts const &gped_covs,
                   ReadStats const &read_stats, Ploidy const ploidy);

    gt_sites const& get_genotyped_records() const {return genotyped_records;}

    static CovCount find_minimum_non_error_cov(double mean_pb_error, poisson_pmf_ptr poisson_prob);
    static likelihood_related_stats make_l_stats(double mean_cov_depth, double mean_pb_error);
    AlleleIds get_haplogroups_with_sites(Marker const& site_ID, AlleleIds candidate_haplogroups) const;
    void invalidate_if_needed(Marker const& parent_site_ID, AlleleIds haplogroups);
};
}

#endif //GRAMTOOLS_INFER_HPP
