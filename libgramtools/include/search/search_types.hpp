/** @file
 * Defines the key data structures supporting `quasimap`ping.
 */
#include "prg/prg.hpp"


#ifndef GRAMTOOLS_SEARCH_TYPES_HPP
#define GRAMTOOLS_SEARCH_TYPES_HPP

// We need a signifier for SearchState with several alleles in the same site
// This signifier must NEVER be a possible allele ID
#define ALLELE_UNKNOWN 0

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
     * Boils down to an SA interval (`gram::SA_Interval`) and a set of variants traversed,
     * currently in traversal and so far (`gram::VariantSitePath`).
     * The former gets used for extending the search while the latter gets used to record coverage information.
     */
    struct SearchState {
        SA_Interval sa_interval = {}; /**< Stores an interval in the suffix array. (By definition,) All members of the interval share a certain prefix of a suffix of the prg.*/
        VariantSitePath traversed_path = {}; /**< Stores the loci that have been entered AND exited during search */
        VariantSitePath traversing_path = {}; /**< Stores the loci that have been entered but not (yet, or ever) exited*/
        SearchVariantSiteState variant_site_state = SearchVariantSiteState::unknown;

        bool invalid = false; /**< Represents whether no path is found in the prg. */

        bool operator==(const SearchState &other) const {
            return this->sa_interval == other.sa_interval
                   and this->traversed_path == other.traversed_path
                   and this->variant_site_state == other.variant_site_state;
        };
    };

    using SearchStates = std::list<SearchState>;
}

#endif //GRAMTOOLS_SEARCH_TYPES_HPP
