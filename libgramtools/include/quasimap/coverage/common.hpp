#include "search/search_types.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_COVERAGE_COMMON_HPP
#define GRAMTOOLS_COVERAGE_COMMON_HPP

namespace coverage {
    namespace record {
        void search_states(Coverage &coverage,
                           const SearchStates &search_states,
                           const uint64_t &read_length,
                           const PRG_Info &prg_info,
                           const uint32_t &random_seed = 0);
    }

    namespace generate {
        Coverage empty_structure(const PRG_Info &prg_info);
    }

    namespace dump {
        void all(const Coverage &coverage,
                 const Parameters &parameters);
    }
}

bool check_allele_encapsulated(const SearchState &search_state,
                               const uint64_t &read_length,
                               const PRG_Info &prg_info);

bool multiple_allele_encapsulated(const SearchState &search_state,
                                  const uint64_t &read_length,
                                  const PRG_Info &prg_info);

#endif //GRAMTOOLS_COVERAGE_COMMON_HPP
