#include <cassert>
#include <fstream>
#include <vector>

#include "search.hpp"

#include "quasimap/utils.hpp"
#include "quasimap/coverage/allele_sum.hpp"


AlleleSumCoverage coverage::generate::allele_sum_structure(const PRG_Info &prg_info) {
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


void coverage::record::allele_sum(Coverage &coverage,
                                  const SearchStates &search_states) {
    auto &allele_sum_coverage = coverage.allele_sum_coverage;

    for (const auto &search_state: search_states) {
        for (const auto &variant_site: search_state.variant_site_path) {
            auto marker = variant_site.first;
            auto allell_id = variant_site.second;

            auto min_boundary_marker = 5;
            auto site_coverage_index = (marker - min_boundary_marker) / 2;
            auto allele_coverage_index = allell_id - 1;

            allele_sum_coverage[site_coverage_index][allele_coverage_index] += 1;
        }
    }
}
