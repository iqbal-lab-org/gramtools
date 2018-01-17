#include <cassert>
#include <fstream>
#include <vector>

#include "search.hpp"

#include "quasimap/utils.hpp"
#include "quasimap/coverage/allele_sum.hpp"


AlleleSumCoverage generate_allele_sum_coverage_structure(const PRG_Info &prg_info) {
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