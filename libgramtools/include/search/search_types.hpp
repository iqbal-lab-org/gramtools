/** @file
 * Defines the key data structures supporting `quasimap`ping.
 */
#include "prg/prg.hpp"


#ifndef GRAMTOOLS_SEARCH_TYPES_HPP
#define GRAMTOOLS_SEARCH_TYPES_HPP

namespace gram {
    /**
     * Expresses the positioning of the current search state relative to variant sites.
     * Initialised at `unknown`.
     */
    enum class SearchVariantSiteState {
        within_variant_site,
        outside_variant_site,
        unknown
    };


    /**
     * A single path of a read through the prg.
     * Boils down to an SA interval (`gram::SA_Interval`) and a set of variants traversed so far (`gram::variant_site_path`).
     * The former gets used for extending the search while the latter gets used to record coverage information.
     */
    struct SearchState {
        SA_Interval sa_interval = {}; /**< Stores an interval in the suffix array. (By definition,) All members of the interval share a certain prefix of a suffix of the prg.*/
        VariantSitePath variant_site_path = {}; /**< Stores the path taken through variant sites in the prg. */
        SearchVariantSiteState variant_site_state = SearchVariantSiteState::unknown;

        bool invalid = false; /**< Represents whether no path is found in the prg. */

        bool operator==(const SearchState &other) const {
            return this->sa_interval == other.sa_interval
                   and this->variant_site_path == other.variant_site_path
                   and this->variant_site_state == other.variant_site_state;
        };
    };

    using SearchStates = std::list<SearchState>;
}

#endif //GRAMTOOLS_SEARCH_TYPES_HPP
