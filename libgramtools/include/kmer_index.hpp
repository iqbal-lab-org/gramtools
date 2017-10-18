#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "utils.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_HPP
#define GRAMTOOLS_KMER_INDEX_HPP

template<typename SEQUENCE>
struct seq_hash {
    std::size_t operator()(const SEQUENCE &seq) const {
        std::size_t hash = 0;
        boost::hash_range(hash, seq.begin(), seq.end());
        return hash;
    }
};

template<typename SEQUENCE, typename T>
using sequence_map = std::unordered_map<SEQUENCE, T, seq_hash<SEQUENCE>>;

template<typename SEQUENCE>
using sequence_set = std::unordered_set<SEQUENCE, seq_hash<SEQUENCE>>;

using SA_Intervals = std::list<SA_Interval>;
using KmerSA_Intervals = sequence_map<Pattern, SA_Intervals>;
using KmerVariantSitePaths = sequence_map<Pattern, VariantSitePaths>;
// Patterns added when not in variant site or when entierly within a single allele
using NonSiteCrossingKmers = sequence_set<Pattern>;

struct KmerIndex {
    KmerSA_Intervals sa_intervals_map;
    KmerVariantSitePaths variant_site_paths_map;
    NonSiteCrossingKmers non_site_crossing_kmers;
};

#endif //GRAMTOOLS_KMER_INDEX_HPP
