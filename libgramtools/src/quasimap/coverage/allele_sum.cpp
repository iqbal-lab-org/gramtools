#include <cassert>
#include <fstream>
#include <vector>

#include "search/search.hpp"

#include "quasimap/utils.hpp"
#include "quasimap/coverage/allele_sum.hpp"


using namespace gram;


AlleleSumCoverage gram::coverage::generate::allele_sum_structure(const PRG_Info &prg_info) {
    uint64_t numer_of_variant_sites = get_number_of_variant_sites(prg_info);
    AlleleSumCoverage allele_sum_coverage(numer_of_variant_sites);

    const auto min_boundary_marker = 5;
    bool last_char_was_zero = true;

    for (const auto &mask_value: prg_info.sites_mask) {
        if (mask_value == 0) {
            last_char_was_zero = true;
            continue;
        }

        const auto &current_marker = mask_value;
        if (last_char_was_zero) {
            auto variant_site_cover_index = (current_marker - min_boundary_marker) / 2;
            allele_sum_coverage[variant_site_cover_index].push_back(0);
            last_char_was_zero = false;
        }
    }
    return allele_sum_coverage;
}


void gram::coverage::record::allele_sum(Coverage &coverage,
                                        const SearchStates &search_states) {
    auto &allele_sum_coverage = coverage.allele_sum_coverage;
    HashSet<VariantSite> seen_sites;

    for (const auto &search_state: search_states) {
        for (const auto &variant_site: search_state.variant_site_path) {
            bool site_seen_previously = seen_sites.find(variant_site) != seen_sites.end();
            if (site_seen_previously)
                continue;

            auto marker = variant_site.first;
            auto allell_id = variant_site.second;

            auto min_boundary_marker = 5;
            auto site_coverage_index = (marker - min_boundary_marker) / 2;
            auto allele_coverage_index = allell_id - 1;

            #pragma omp atomic
            allele_sum_coverage[site_coverage_index][allele_coverage_index] += 1;
            seen_sites.insert(variant_site);
        }
    }
}


void gram::coverage::dump::allele_sum(const Coverage &coverage,
                                      const Parameters &parameters) {
    std::ofstream file_handle(parameters.allele_sum_coverage_fpath);
    for (const auto &variant_site_coverage: coverage.allele_sum_coverage) {
        auto allele_count = 0;
        for (const auto &sum_coverage: variant_site_coverage) {
            file_handle << sum_coverage;
            auto not_last_coverage = allele_count++ < variant_site_coverage.size() - 1;
            if (not_last_coverage)
                file_handle << " ";
        }
        file_handle << std::endl;
    }
}