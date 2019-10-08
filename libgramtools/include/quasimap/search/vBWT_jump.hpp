#ifndef GRAMTOOLS_vBWT_JUMP
#define GRAMTOOLS_vBWT_JUMP

#include "common/utils.hpp"
#include "quasimap/search_types.hpp"
#include <vector>

namespace gram {
    using MarkersSearchResults = std::vector<VariantLocus>;
    using Locus_and_SearchState = std::pair<VariantLocus, SearchState>;
    using Locus_and_SearchStates = std::vector<Locus_and_SearchState>;

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

    /**
     * This function finds all variant markers (site or allele) inside the BWT within a given SA interval.
     * Indeed, if a variant marker precedes an index position of the SA interval (as discovered using BWT),
     * the search states will need to be updated accordingly.
     *
     * @return A vector of `VariantLocus`
     */
    MarkersSearchResults left_markers_search(const SearchState &search_state,
                                             const PRG_Info &prg_info);

    Locus_and_SearchState extend_targets_site_exit(VariantLocus const &target_locus, SearchState const &search_state,
                                                   PRG_Info const &prg_info);

    Locus_and_SearchStates extend_targets_site_entry(VariantLocus const &target_locus, SearchState const &search_state,
                                                     PRG_Info const &prg_info);

/**
 * Computes the full SA interval of a given allele marker.
 * Note that the way this is computed is robust to variant markers not being continuous:
 *      eg, could have site with markers 5/6 & another with 9/10 without a site with 7/8.
 *
 * `sigma`: the alphabet size
 * `char2comp`: maps a symbol of the alphabet to a positive integer in [0...sigma - 1]
 * `comp2char`: the inverse mapping of above
 */
    SA_Interval get_allele_marker_sa_interval(const Marker &allele_marker_char,
                                              const PRG_Info &prg_info);
}

#endif //GRAMTOOLS_vBWT_JUMP
