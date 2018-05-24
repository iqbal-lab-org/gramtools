#include "search/search_types.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_ALLELE_BASE_HPP
#define GRAMTOOLS_ALLELE_BASE_HPP

namespace coverage {
    namespace generate {
        SitesAlleleBaseCoverage allele_base_structure(const PRG_Info &prg_info);
    }

    namespace record {
        void allele_base(Coverage &coverage,
                         const SearchStates &search_states,
                         const uint64_t &read_length,
                         const PRG_Info &prg_info);
    }

    namespace dump {
        void allele_base(const Coverage &coverage,
                         const Parameters &parameters);
    }
}

std::string dump_allele_base_coverage(const SitesAlleleBaseCoverage &sites);

uint64_t inter_site_base_count(const uint64_t &first_site_marker,
                               const uint64_t &second_site_marker,
                               const PRG_Info &prg_info);

using SitesCoverageBoundaries = PairHashMap<VariantSite, uint64_t>;

uint64_t set_site_base_coverage(Coverage &coverage,
                                SitesCoverageBoundaries &sites_coverage_boundaries,
                                const VariantSite &path_element,
                                const uint64_t allele_coverage_offset,
                                const uint64_t max_bases_to_set);

uint64_t allele_start_offset_index(const uint64_t within_allele_prg_index, const PRG_Info &prg_info);

#endif //GRAMTOOLS_ALLELE_BASE_HPP
