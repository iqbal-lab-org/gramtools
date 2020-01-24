#include <cmath>
#include "genotype/infer/infer.hpp"
#include "genotype/infer/allele_extracter.hpp"

likelihood_related_stats
LevelGenotyper::make_l_stats(double mean_cov_depth, double mean_pb_error){
    PoissonLogPmf poisson_prob{params{mean_cov_depth}};

    // store natural log of pb error also because of its use in likelihood formulae
    return likelihood_related_stats{
        mean_cov_depth,
        mean_pb_error, log(mean_pb_error),
        log(1 - exp(-mean_cov_depth)),
        log(1 - exp(-mean_cov_depth / 2)),
        find_minimum_non_error_cov(mean_pb_error, &poisson_prob),
        poisson_prob,
        PoissonLogPmf{params{mean_cov_depth / 2}}
    };
}

CovCount LevelGenotyper::find_minimum_non_error_cov(double mean_pb_error, poisson_pmf_ptr poisson_prob) {
    CovCount min_count{1};
    while ( (*poisson_prob)(params{static_cast<double>(min_count)}) <= log(pow(mean_pb_error, min_count)) )
        ++min_count;
    return min_count;
}

LevelGenotyper::LevelGenotyper(coverage_Graph const& cov_graph, SitesGroupedAlleleCounts const& gped_covs,
                               ReadStats const& read_stats) :
        cov_graph(&cov_graph), gped_covs(&gped_covs){
    genotyped_records.reserve(cov_graph.bubble_map.size()); // Pre-allocate one slot for each bubble in the PRG
    auto mean_cov_depth = read_stats.get_mean_cov_depth();
    auto mean_pb_error = read_stats.get_mean_pb_error();

    PoissonLogPmf poisson_prob{params{mean_cov_depth}};

    l_stats = std::move(make_l_stats(mean_cov_depth, mean_pb_error));

    // Genotype each bubble in the PRG, in most nested to less nested order.
    for (auto const& bubble_pair : cov_graph.bubble_map){
        auto site_index = siteID_to_index(bubble_pair.first->get_site_ID());
        // TODO : set end node of site here
        auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
        auto extracted_alleles = extracter.get_alleles();
        auto& gped_covs_for_site = gped_covs.at(site_index);
        auto genotyped = LevelGenotyperModel(&extracted_alleles, &gped_covs_for_site,
                Ploidy::Haploid, &l_stats);

        genotyped_records.at(site_index) = genotyped.get_site();

        // TODO: invalidation
    }
}
