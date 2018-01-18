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
}

#endif //GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP
