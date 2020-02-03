#include <cmath>
#include "genotype/infer/level_genotyping/runner.hpp"
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

AlleleIds LevelGenotyper::get_haplogroups_with_sites(Marker const& site_ID, AlleleIds candidate_haplogroups) const{
    AlleleIds result{};
    if (child_m.find(site_ID) == child_m.end()) return result;
    auto child_entry = child_m.at(site_ID);
    for (auto const& candidate : candidate_haplogroups){
       if (child_entry.find(candidate) != child_entry.end()) result.push_back(candidate);
    }
    return result;
}

void LevelGenotyper::invalidate_if_needed(Marker const& parent_site_ID, AlleleIds haplogroups){
    if (haplogroups.size() == 0) return;

    std::vector<VariantLocus> to_process;
    for (auto const& haplogroup : haplogroups) to_process.push_back(VariantLocus{parent_site_ID, haplogroup});

    VariantLocus cur_locus;
    while (to_process.size() > 0){
        cur_locus = to_process.back();
        to_process.pop_back();

        // We know the site/haplogroup combination bears 1+ site so use .at() .
        auto sites_on_haplogroup = child_m.at(cur_locus.first).at(cur_locus.second);
        for (auto const& child_marker : sites_on_haplogroup){
            auto referent_site = genotyped_records.at(siteID_to_index(child_marker));
            if (referent_site->is_null()) continue;
            referent_site->make_null();

            auto all_haplos = referent_site->get_all_haplogroups();
            auto haplos_with_sites = get_haplogroups_with_sites(child_marker, all_haplos);
            for (auto const& h : haplos_with_sites) to_process.push_back(VariantLocus{child_marker,h});
        }
    }
}

LevelGenotyper::LevelGenotyper(coverage_Graph const &cov_graph, SitesGroupedAlleleCounts const &gped_covs,
                               ReadStats const &read_stats, Ploidy const ploidy) :
        cov_graph(&cov_graph), gped_covs(&gped_covs), ploidy(ploidy){
    child_m = build_child_map(cov_graph.par_map); // Required for site invalidation
    genotyped_records.resize(cov_graph.bubble_map.size()); // Pre-allocate one slot for each bubble in the PRG

    auto mean_cov_depth = read_stats.get_mean_cov_depth();
    auto mean_pb_error = read_stats.get_mean_pb_error();
    PoissonLogPmf poisson_prob{params{mean_cov_depth}};
    l_stats = std::move(make_l_stats(mean_cov_depth, mean_pb_error));

    // Genotype each bubble in the PRG, in most nested to less nested order.
    for (auto const& bubble_pair : cov_graph.bubble_map){
        auto site_ID = bubble_pair.first->get_site_ID();
        auto site_index = siteID_to_index(site_ID);

        auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
        auto& gped_covs_for_site = gped_covs.at(site_index);
        auto extracted_alleles = extracter.get_alleles();

        auto genotyped = LevelGenotyperModel(&extracted_alleles, &gped_covs_for_site,
                                             ploidy, &l_stats, ! extracter.ref_allele_got_made_naturally());
        auto genotyped_site = genotyped.get_site();
        genotyped_records.at(site_index) = genotyped_site;
        // Line below is so that when allele extraction occurs and jumps through a previously
        // genotyped site, it knows where in the graph to resume from.
        genotyped_records.at(site_index)->set_site_end_node(bubble_pair.second);

        // Invalidation process attempted only if this site contains 1+ site
        if (! genotyped_site->is_null() && child_m.find(site_ID) != child_m.end()){
            auto candidate_haplogroups =
                    genotyped_site->get_nonGenotyped_haplogroups();
            auto haplogroups_with_sites = get_haplogroups_with_sites(site_ID, candidate_haplogroups);
            invalidate_if_needed(site_ID, haplogroups_with_sites);
        }
    }
}
