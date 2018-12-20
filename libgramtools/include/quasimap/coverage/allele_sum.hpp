/** @file
 * Defines coverage related operations for allele sum coverage.
 * `AlleleSumCoverage` stores the sum of all reads mapped for each allele of each variant site.
 */
#include "search/search_types.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_ALLELE_SUM_HPP
#define GRAMTOOLS_ALLELE_SUM_HPP

namespace gram::coverage {
    namespace generate {
        AlleleSumCoverage allele_sum_structure(const PRG_Info &prg_info);
    }

    namespace record {
        /**
         * Loops over search_states and variant site paths, recording traversed alleles.
         * @param coverage The `Coverage` structure common to all mapped reads.
         * @param search_states The selected `SearchStates` for recording coverage.
         */
        void allele_sum(Coverage &coverage,
                        const SearchStates &search_states);
    }

    namespace dump {
        void allele_sum(const Coverage &coverage,
                        const Parameters &parameters);
    }
}

#endif //GRAMTOOLS_ALLELE_SUM_HPP
