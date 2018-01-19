#include <vector>
#include <list>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "sequence_read/seqread.hpp"


#ifndef GRAMTOOLS_UTILS_HPP
#define GRAMTOOLS_UTILS_HPP

template<typename SEQUENCE>
struct seq_hash {
    std::size_t operator()(const SEQUENCE &seq) const {
        std::size_t hash = 0;
        boost::hash_range(hash, seq.begin(), seq.end());
        return hash;
    }
};

template<typename SEQUENCE, typename T>
using SequenceHashMap = std::unordered_map<SEQUENCE, T, seq_hash<SEQUENCE>>;

template<typename T>
using HashSet = std::unordered_set<T, boost::hash<T>>;

using Base = uint8_t;
using Pattern = std::vector<Base>;
using Patterns = std::vector<Pattern>;

using Marker = uint64_t;
using AlleleId = uint64_t;

using VariantSite = std::pair<Marker, AlleleId>;
using VariantSitePath = std::list<VariantSite>;
using VariantSitePaths = std::list<VariantSitePath>;

using SA_Index = uint64_t;
using SA_Interval = std::pair<SA_Index, SA_Index>;

Pattern reverse_compliment_read(const Pattern &read);
Pattern encode_dna_bases(const std::string &dna_str);
Pattern encode_dna_bases(const GenomicRead &read_sequence);

Base encode_dna_base(const char &base_str);

#endif //GRAMTOOLS_UTILS_HPP
