#include <fstream>
#include <vector>

#include "search.hpp"

#include "quasimap/utils.hpp"
#include "quasimap/coverage/grouped_allele_counts.hpp"


SitesGroupedAlleleCounts coverage::generate::grouped_allele_counts(const PRG_Info &prg_info) {
    uint64_t numer_of_variant_sites = get_number_of_variant_sites(prg_info);
    SitesGroupedAlleleCounts grouped_allele_counts(numer_of_variant_sites);
    return grouped_allele_counts;
}


void coverage::record::grouped_allele_counts(Coverage &coverage,
                                             const SearchStates &search_states) {
    using AlleleIdSet = std::set<AlleleId>;
    std::unordered_map<Marker, AlleleIdSet> site_allele_group;

    for (const auto &search_state: search_states) {
        for (const auto &variant_site: search_state.variant_site_path) {
            auto site_marker = variant_site.first;
            auto allell_id = variant_site.second;
            site_allele_group[site_marker].insert(allell_id);
        }
    }

    for (const auto &entry: site_allele_group) {
        auto site_marker = entry.first;
        auto allele_ids_set = entry.second;
        AlleleIds allele_ids;
        std::copy(allele_ids_set.begin(), allele_ids_set.end(),
                  std::back_inserter(allele_ids));

        auto min_boundary_marker = 5;
        auto site_coverage_index = (site_marker - min_boundary_marker) / 2;

        auto &site_coverage = coverage.grouped_allele_counts[site_coverage_index];
        site_coverage[allele_ids] += 1;
    }
}