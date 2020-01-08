/** @file
 * Defines the kmer index and the caching structure for remembering the relevant previous mappings.
 */
#include "common/utils.hpp"
#include "genotype/quasimap/search_types.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_TYPES_HPP
#define GRAMTOOLS_KMER_INDEX_TYPES_HPP

// s1 s2 s3 s4
// 1  2   1  4

namespace gram {
    struct CacheElement {
        SearchStates search_states = {};
        int_Base base = 0; /**< The next base added relative to the previous `CacheElement`.*/
    };
    using KmerIndexCache = std::list<CacheElement>; /**< Stored previously computed `SearchStates` for re-use when indexing different kmers. */

    using KmerIndex = SequenceHashMap<Pattern, SearchStates>; /**< Unordered map linking a kmer to all mapped locations in the prg.*/
}

#endif //GRAMTOOLS_KMER_INDEX_TYPES_HPP
