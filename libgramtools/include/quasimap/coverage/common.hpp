#include "search.hpp"
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

    void dump(const Coverage &coverage,
              const Parameters &parameters);
}

SA_Interval random_select_sa_interval(const SearchStates &search_states,
                                      const uint32_t &random_seed = 0);

SearchStates filter_for_sa_interval(const SA_Interval &target_sa_interval,
                                    const SearchStates &search_states);

#endif //GRAMTOOLS_COVERAGE_COMMON_HPP
