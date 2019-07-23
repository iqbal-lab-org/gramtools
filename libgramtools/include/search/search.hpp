/** @file
 * Procedures supporting variant aware backward searching through the prg.
 * @note `char2comp` attribute of `fm_index` gives the lexicographic ordering of the queried symbol. This allows for
 * finding symbol's first occurrence in the SA using the `C` array. For eg, we do not assume that site marker '5' is
 * the 5th element of the `C` array, because we can be given a prg which can have discontinuous integers marking variant sites.
 */
#include "common/utils.hpp"
#include "kmer_index/kmer_index_types.hpp"
#include "search_types.hpp"


#ifndef GRAMTOOLS_SEARCH_HPP
#define GRAMTOOLS_SEARCH_HPP

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
     * Generates a list of `SearchState`s from a read and a kmer, which is 3'-most kmer in the read.
     * The kmer_index is queried to generate an initial set of `SearchState`s (precomputed at `build` stage) to start from.
     * @return SearchStates: a list of `SearchState`s, which at core are an SA interval and a path through the prg (marker-allele ID pairs)
     */
    SearchStates search_read_backwards(const Pattern &read,
                                       const Pattern &kmer,
                                       const KmerIndex &kmer_index,
                                       const PRG_Info &prg_info);

    /**
     * Updates each SearchState with the next character in the read.
     * @param pattern_char the next character in the read to look for in the prg.
     * @param SearchStates a set of SearchState elements; each contains an SA interval.
     */
    SearchStates search_base_backwards(const Base &pattern_char,
                                       const SearchStates &search_states,
                                       const PRG_Info &prg_info);
using SaIndexRightOfMarker = uint64_t;
    using MarkersSearchResult = std::pair<SaIndexRightOfMarker, Marker>; /** SaIndexRightOfMarker is the SA index position of the character just right of the marker in the prg.*/
    using MarkersSearchResults = std::vector<MarkersSearchResult>;

    /**
     * This function finds all variant markers (site or allele) inside the BWT within a given SA interval.
     * Indeed, if a variant marker precedes an index position of the SA interval (as discovered using BWT),
     * the search states will need to be updated accordingly.
     *
     * @return A vector of markers, which are pairs of SA index position and marker character.
     * @note The SA index position is not that of the marker, but of the character to the right of it in the prg.
     */
    MarkersSearchResults left_markers_search(const SearchState &search_state,
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

    /**
     * On leaving a site, search states have cached path elements. These path elements describe where the search state has
     * just been. We handle/flush the cache once we've dealt with a character outside of the site or the last character of
     * the read is encountered.
     *
     * The cache consists of a single path element (site marker, allele id).
     */
    void process_search_state_path_cache(SearchState &search_state);

    /**
     * **The key read mapping procedure**.
     * First updates SA_intervals to search next based on variant marker presence.
     * Then executes regular backward search.
     */
    SearchStates process_read_char_search_states(const Base &pattern_char,
                                                 const SearchStates &old_search_states,
                                                 const PRG_Info &prg_info);

    /**
    * Retrieve the allele id from an allele marker SA index.
    * The previous position in the original text is used to query the `allele_mask`, part of `PRG_Info`. It contains
    * the allele id for each position of the prg.
    */
    AlleleId get_allele_id(const SA_Index &allele_marker_sa_index,
                           const PRG_Info &prg_info);

    /**
    * Finds the whole SA interval associated with a given allele marker.
    * Remember: they are contiguous in the SA.
    * The procedure queries the C array for the first occurrence of the marker allele, and the first occurrence of the next
    * site marker.
    */
    SA_Interval get_allele_marker_sa_interval(const Marker &site_marker_char,
                                              const PRG_Info &prg_info);

    /**
    * Call `process_markers_search_state` for each `SearchState`.
    * Each SA index whose corresponding BWT entry is a marker will generate one or more new `SearchStates`.
    * Note that the original `SearchState` is otherwise left untouched; SA indices with preceding markers in the prg will get naturally dropped by backward base extension.
    * @see process_markers_search_state()
    * @see SearchState()
    */
    SearchStates process_markers_search_states(const SearchStates &search_states,
                                               const PRG_Info &prg_info);

    /**
     * For a given `SearchState`, add new `SearchState`s based on variant marker presence.
     * Variant markers are found by querying the BWT on the SA interval of the `SearchState`.
     * New `SearchState`s are then generated based on whether site or allele markers are found.
     * @see left_markers_search()
     *
     */
    SearchStates process_markers_search_state(const SearchState &search_state,
                                              const PRG_Info &prg_info);

    std::string serialize_search_state(const SearchState &search_state);

    std::ostream &operator<<(std::ostream &os, const SearchState &search_state);

}

#endif //GRAMTOOLS_SEARCH_HPP
