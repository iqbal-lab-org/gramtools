#include "search.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP
#define GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP

namespace coverage {
    namespace generate {
        SitesGroupedAlleleCounts grouped_allele_counts(const PRG_Info &prg_info);
    }

    namespace record {
        void grouped_allele_counts(Coverage &coverage,
                                   const SearchStates &search_states);
    }

    namespace dump {
        void grouped_allele_counts(const Coverage &coverage,
                                   const Parameters &parameters);
    }
}

AlleleGroupHash hash_allele_groups(const SitesGroupedAlleleCounts &sites);

std::string dump_site(const AlleleGroupHash &allele_ids_groups_hash,
                      const GroupedAlleleCounts &site);

std::string dump_site_counts(const AlleleGroupHash &allele_ids_groups_hash,
                             const SitesGroupedAlleleCounts &sites);

std::string dump_allele_groups(const AlleleGroupHash &allele_ids_groups_hash);

std::string dump_grouped_allele_counts(const SitesGroupedAlleleCounts &sites);

#endif //GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP
