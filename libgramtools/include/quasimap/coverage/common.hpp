#include "quasimap/search_types.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_COVERAGE_COMMON_HPP
#define GRAMTOOLS_COVERAGE_COMMON_HPP

namespace gram {

    /**
     * Each type of coverage operation (record, generate, dump) operates on each level of coverage information.
     */
    namespace coverage {
        namespace record {
            /**
             * Selects read mappings and records coverage information.
             * @see selection()
             */
            void search_states(Coverage &coverage,
                               const SearchStates &search_states,
                               const uint64_t &read_length,
                               const PRG_Info &prg_info,
                               const uint32_t &random_seed = 0);
        }

        namespace generate {
            /**
             * Calls the routines for building empty structures to record different types of coverage information.
             */
            Coverage empty_structure(const PRG_Info &prg_info);
        }

        namespace dump {
            /**
             * Write coverage information to disk.
             */
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

    uint64_t random_int_inclusive(const uint64_t &min,
                                  const uint64_t &max,
                                  const uint64_t &random_seed);

    uint32_t count_nonvariant_search_states(const SearchStates &search_states);

    using SitePath = std::set<Marker>;
    using uniqueSitePaths = std::map<SitePath, SearchStates>;

    SitePath get_path_sites(const SearchState &search_state);

    uniqueSitePaths get_unique_site_paths(const SearchStates &search_states);

}

#endif //GRAMTOOLS_COVERAGE_COMMON_HPP
