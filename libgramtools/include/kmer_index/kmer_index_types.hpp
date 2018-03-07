#include "common/utils.hpp"
#include "search/search_types.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_TYPES_HPP
#define GRAMTOOLS_KMER_INDEX_TYPES_HPP

struct CacheElement {
    SearchStates search_states = {};
    Base base = 0;
};
using KmerIndexCache = std::list<CacheElement>;

using KmerIndex = SequenceHashMap<Pattern, SearchStates>;

#endif //GRAMTOOLS_KMER_INDEX_TYPES_HPP
