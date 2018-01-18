#include "search.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_ALLELE_SUM_HPP
#define GRAMTOOLS_ALLELE_SUM_HPP

namespace coverage {
    namespace generate {
        AlleleSumCoverage allele_sum_structure(const PRG_Info &prg_info);
    }

    namespace record {
        void allele_sum(Coverage &coverage,
                        const SearchStates &search_states);
    }
}

#endif //GRAMTOOLS_ALLELE_SUM_HPP
