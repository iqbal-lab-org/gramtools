#ifndef GRAMTOOLS_vBWT_JUMP
#define GRAMTOOLS_vBWT_JUMP

#include "common/utils.hpp"
#include "genotype/quasimap/search/types.hpp"
#include <vector>

namespace gram {
    using MarkersSearchResults = std::vector<VariantLocus>;
    /**
     * A convenience structure for holding together a `VariantLocus` that needs a vBWT jump and a
     * SearchState for holding the already traversed path of `VariantLocus`
     */
    struct Locus_and_SearchState{
        VariantLocus locus;
        SearchState search_state;
        bool commit_me;
    };
    using Locus_and_SearchStates = std::vector<Locus_and_SearchState>;

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
     * This function finds all variant markers (site or allele) inside the BWT within a given SA interval.
     * Indeed, if a variant marker precedes an index position of the SA interval,
     * the search states will need to be updated accordingly.
     *
     * @return A vector of `VariantLocus`
     */
    MarkersSearchResults left_markers_search(const SearchState &search_state,
                                             const PRG_Info &prg_info);

    /**
     * For a given `SearchState`, add new `SearchState`s if there are variant markers preceding any position in the
     * SA interval.
     *
     * The function deals with adjacent variant markers and delegates their processing to:
     *      - `extend_targets_site_exit` for leaving a site
     *      - `extend_targets_site_entry` for entering a site
     *
     * Two questions are systematically asked for each new `SearchState`:
     *      1) Does it need to be processed further, due to adjacent variant markers?
     *      2) Does it need to be backward searched?
     */
    SearchStates search_state_vBWT_jumps(const SearchState &current_search_state,
                                         const PRG_Info &prg_info);


    /**
     * We are leaving a site.
     * The contract is as follows:
     *      - This function will return a single locus and search state, because a site entry is a single point
     *      - The locus can be:
     *          - 'empty': the search state will have exited the left-most site marker in case of multiple exits.
     *          - an allele marker: when the exit point is followed by an entry point (allele marker)
     */
    Locus_and_SearchState extend_targets_site_exit(VariantLocus const &target_locus, SearchState const &search_state,
                                                   PRG_Info const &prg_info);

    /**
     * We are entering a site.
     * We create a new `SearchState` able to map into the alleles of the site,
     * and register any adjacent variant markers for further processing (but not backward searching)
     */
    Locus_and_SearchStates extend_targets_site_entry(VariantLocus const &target_locus, SearchState const &search_state,
                                                     PRG_Info const &prg_info);

}

#endif //GRAMTOOLS_vBWT_JUMP
