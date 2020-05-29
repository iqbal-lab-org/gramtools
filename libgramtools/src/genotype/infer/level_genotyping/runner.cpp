#include <cmath>

#include "genotype/infer/level_genotyping/runner.hpp"
#include "genotype/infer/allele_extracter.hpp"
#include "prg/coverage_graph.hpp"
#include "prg/make_data_structures.hpp"
#include "genotype/infer/output_specs/fields.hpp"

using namespace gram::genotype::output_spec;

header_vec LevelGenotyper::get_model_specific_headers(){
    header_vec result{
        vcf_meta_info_line{"Model", "LevelGenotyping"},
        vcf_meta_info_line{"FORMAT", "GT_CONF",
                           "Genotype confidence as "
                 "likelihood ratio of called and next most likely genotype.",
                 "1", "Float"}
    };
    return result;
}

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



LevelGenotyper::LevelGenotyper(coverage_Graph const &cov_graph, SitesGroupedAlleleCounts const &gped_covs,
                               ReadStats const &read_stats, Ploidy const ploidy) :
                               ploidy(ploidy){
    this->cov_graph = &cov_graph;
    this->gped_covs = &gped_covs;
    child_m = build_child_map(cov_graph.par_map); // Required for site invalidation
    genotyped_records.resize(cov_graph.bubble_map.size()); // Pre-allocate one slot for each bubble in the PRG

    auto mean_cov_depth = read_stats.get_mean_cov_depth();
    auto mean_pb_error = read_stats.get_mean_pb_error();
    l_stats = std::move(make_l_stats(mean_cov_depth, mean_pb_error));

    // Genotype each bubble in the PRG, in most nested to less nested order.
    for (auto const& bubble_pair : cov_graph.bubble_map){
        auto site_ID = bubble_pair.first->get_site_ID();
        auto site_index = siteID_to_index(site_ID);

        auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
        auto extracted_alleles = extracter.get_alleles();
        auto& gped_covs_for_site = gped_covs.at(site_index);

        ModelData data(&extracted_alleles, &gped_covs_for_site,
                       ploidy, &l_stats, ! extracter.ref_allele_got_made_naturally());
        auto genotyped = LevelGenotyperModel(data);
        auto genotyped_site = genotyped.get_site();
        genotyped_site->set_pos(bubble_pair.first->get_pos());
        genotyped_records.at(site_index) = genotyped_site;
        // Line below is so that when allele extraction occurs and jumps through a previously
        // genotyped site, it knows where in the graph to resume from.
        genotyped_records.at(site_index)->set_site_end_node(bubble_pair.second);

        run_invalidation_process(genotyped_site, site_ID);
    }
}

