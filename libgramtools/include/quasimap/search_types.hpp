/** @file
 * Defines the key data structures supporting `quasimap`ping.
 */
#include "prg/prg_info.hpp"


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

    /** A path through variant sites is a list of allele/site combinations. */
    using VariantSitePath = std::vector<VariantLocus>;
    using VariantSitePaths = std::vector<VariantSitePath>;

    /** The suffix array (SA) holds the starting index of all (lexicographically sorted) cyclic permutations of the prg.
     * An `SA_Index` is an index into one such position.*/
    using SA_Index = uint32_t;
    using SA_Interval = std::pair<SA_Index, SA_Index>; /**< A set of **contiguous** indices in the suffix array.*/

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

        /**
         * Asks if `SearchState` has crossed any site boundary markers.
         * If it has not it may still have mapped fully inside an allele.
         */
        bool has_path() const{
            return (! this->traversed_path.empty() || ! this->traversing_path.empty());
        }
    };

    using SearchStates = std::list<SearchState>;
}

#endif //GRAMTOOLS_SEARCH_TYPES_HPP
