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

#define CONF_DISTRIB_SIZE 10000

#include "site.hpp"
#include "genotype/parameters.hpp"
#include "probabilities.hpp"
#include "genotype/read_stats.hpp"

using namespace gram::genotype::infer;
using namespace gram::genotype::infer::probabilities;

namespace gram::genotype::infer {

using lvlgt_site_ptr = std::shared_ptr<LevelGenotypedSite>;

class LevelGenotyper : public Genotyper {
    likelihood_related_stats l_stats;
    Ploidy ploidy;

public:
    LevelGenotyper() = default;
    LevelGenotyper(child_map const& ch, gt_sites const& sites) : Genotyper(sites, ch) {}

    LevelGenotyper(coverage_Graph const &cov_graph, SitesGroupedAlleleCounts const &gped_covs,
                   ReadStats const &read_stats, Ploidy ploidy, bool get_gcp = false);

    header_vec get_model_specific_headers() override;

    /**
     * Invalidation of non-chosen haplogroups in nested prgs
     */
    AlleleIds get_haplogroups_with_sites(Marker const& site_ID, AlleleIds candidate_haplogroups) const;
    void invalidate_if_needed(Marker const& parent_site_ID, AlleleIds haplogroups);
    void run_invalidation_process(lvlgt_site_ptr const& genotyped_site, Marker const& site_ID);

    std::vector<double> static get_gtconf_distrib(gt_sites const& input_sites,
                                                  likelihood_related_stats const& input_lstats,
                                                  Ploidy const& input_ploidy);
    static CovCount find_minimum_non_error_cov(double mean_pb_error, pmf_ptr pmf);
    static likelihood_related_stats make_l_stats(double mean_cov, double var_cov, double mean_pb_error);

};
}

#endif //LVLGT_RUNNER
