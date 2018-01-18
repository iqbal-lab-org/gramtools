#include "search.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_COVERAGE_COMMON_HPP
#define GRAMTOOLS_COVERAGE_COMMON_HPP

namespace coverage {
    namespace record {
        void search_states(Coverage &coverage,
                           const SearchStates &search_states,
                           const uint64_t &read_length,
                           const PRG_Info &prg_info);
    }
    
    namespace generate {
        Coverage empty_structure(const PRG_Info &prg_info);
    }

    void dump(const Coverage &coverage,
              const Parameters &parameters);
}

#endif //GRAMTOOLS_COVERAGE_COMMON_HPP
