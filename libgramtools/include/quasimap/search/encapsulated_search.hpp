#ifndef ENCAPS_SEARCH_HPP
#define ENCAPS_SEARCH_HPP

#include "quasimap/search_types.hpp"
namespace gram {

    /**
     * Potentially splits a search state based on whether it is encapsulated within an allele.
     * By doing so, we can assign paths to search states which were previously unknown.
     *
     * Furthermore, it splits search states which are outside of alleles based on SA index. This ensures that the total
     * number of search states allows for the deliberate assignment of coverage (random sampling of one read for multi-mapped reads).
     *
     * Mappings which are encapsulated and occupy the same allele are represented by a single search state.
     */
    SearchStates handle_allele_encapsulated_state(const SearchState &search_state,
                                                  const PRG_Info &prg_info);

    /**
     * @see handle_allele_encapsulated_state()
     */
    SearchStates handle_allele_encapsulated_states(const SearchStates &search_states,
                                                   const PRG_Info &prg_info);
}

#endif //GRAMTOOLS_SEARCH_HPP
