/** @file
 * Procedures supporting variant aware backward searching through the prg.
 * @note `char2comp` attribute of `fm_index` gives the lexicographic ordering of the queried symbol. This allows for
 * finding symbol's first occurrence in the SA using the `C` array. For eg, we do not assume that site marker '5' is
 * the 5th element of the `C` array, because we can be given a prg which can have discontinuous integers marking variant sites.
 */

#ifndef GRAMTOOLS_SEARCH_HPP
#define GRAMTOOLS_SEARCH_HPP

#include "common/utils.hpp"
#include "kmer_index/kmer_index_types.hpp"
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

    /**
     * Situation: we have fully mapped a read to the PRG.
     * Some `SearchState`s may still have unknown allele ids. Here we set those.
     * Modifies the `SearchStates` in place.
     */
        void set_allele_ids(SearchStates &search_states,
                            const PRG_Info &prg_info);
        
    /**
     * Updates each SearchState with the next character in the read.
     * @param pattern_char the next character in the read to look for in the prg.
     * @param SearchStates a set of SearchState elements; each contains an SA interval.
     */
    SearchStates search_base_backwards(const int_Base &pattern_char,
                                       const SearchStates &search_states,
                                       const PRG_Info &prg_info);

    /**
     * Update the current SA interval to include the next character.
     * This is a backward search. SA interval is updated using rank queries on the bwt.
     * @param next_char the next character to look for.
     * @param next_char_first_sa_index the position of the first occurrence of `next_char` in the SA.
     */
    SA_Interval base_next_sa_interval(const Marker &next_char,
                                      const SA_Index &next_char_first_sa_index,
                                      const SA_Interval &current_sa_interval,
                                      const PRG_Info &prg_info);


    std::string serialize_search_state(const SearchState &search_state);

    std::ostream &operator<<(std::ostream &os, const SearchState &search_state);

}

#endif //GRAMTOOLS_SEARCH_HPP
