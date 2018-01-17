#include <unordered_set>
#include <unordered_map>

#include "utils.hpp"
#include "search_types.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_TYPES_HPP
#define GRAMTOOLS_KMER_INDEX_TYPES_HPP

struct CacheElement {
    SearchStates search_states;
    Base base = 0;
};
using KmerIndexCache = std::list<CacheElement>;

template<typename SEQUENCE, typename T>
using sequence_map = std::unordered_map<SEQUENCE, T, seq_hash<SEQUENCE>>;
using KmerIndex = sequence_map<Pattern, SearchStates>;

#endif //GRAMTOOLS_KMER_INDEX_TYPES_HPP
