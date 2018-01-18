#include <fstream>
#include <vector>

#include "search.hpp"

#include "quasimap/utils.hpp"
#include "quasimap/coverage/grouped_allele_counts.hpp"


SitesGroupedAlleleCounts coverage::generate::grouped_allele_counts(const PRG_Info &prg_info) {
    uint64_t numer_of_variant_sites = get_number_of_variant_sites(prg_info);
    SitesGroupedAlleleCounts grouped_allele_counts;
    grouped_allele_counts.reserve(numer_of_variant_sites);
    return grouped_allele_counts;
}


void coverage::record::grouped_allele_counts(Coverage &coverage,
                                             const SearchStates &search_states) {
    auto &allele_sum_coverage = coverage.allele_sum_coverage;
    for (const auto &search_state: search_states) {
        /*
        for (const auto &variant_site: search_state.variant_site_path) {
            auto marker = variant_site.first;
            auto allell_id = variant_site.second;

            auto min_boundary_marker = 5;
            auto variant_site_coverage_index = (marker - min_boundary_marker) / 2;
            auto allele_coverage_index = allell_id - 1;

            allele_sum_coverage[variant_site_coverage_index][allele_coverage_index] += 1;
        }
        */
    }
}