/** @file
 * Defines a set of data structures and functions common to the backend.
 */
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

namespace gram {
    template<typename SEQUENCE>
    struct sequence_hash {
        std::size_t operator()(const SEQUENCE &seq) const {
            std::size_t hash = 0; // Used as an initial seed
            boost::hash_range(hash, seq.begin(), seq.end());
            return hash;
        }
    };

    template<typename SEQUENCE, typename T>
    using SequenceHashMap = std::unordered_map<SEQUENCE, T, sequence_hash<SEQUENCE>>;

    template<typename PAIR, typename T>
    using PairHashMap = std::unordered_map<PAIR, T, boost::hash<PAIR>>;

    template<typename T>
    using HashSet = std::unordered_set<T, boost::hash<T>>;

    using Base = uint8_t; /**< nucleotide represented as byte-sized integer */
    using Pattern = std::vector<Base>; /** A string of nucleotides is represented as a vector of `Base`s. */
    using Patterns = std::vector<Pattern>;

    using Marker = uint64_t; /**< An integer >=5 describing a site or allele marker in the prg. */
    using AlleleId = uint64_t; /**< An integer describing which allele is referred to within a given variant site. */

    using VariantLocus = std::pair<Marker, AlleleId>; /**< A Variant site/`AlleleId` combination.*/
    using VariantSitePath = std::list<VariantLocus>; /**< A path through variant sites is a list of allele/site combinations. */
    using VariantSitePaths = std::list<VariantSitePath>;

    /** The suffix array (SA) holds the starting index of all (lexicographically sorted) cyclic permutations of the prg.
     * An `SA_Index` is an index into one such position.*/
    using SA_Index = uint64_t;
    using SA_Interval = std::pair<SA_Index, SA_Index>; /**< A set of **contiguous** indices in the suffix array.*/

    /**
     * Produce the reverse complement of a `read`.
     */
    Pattern reverse_complement_read(const Pattern &read);

    Pattern encode_dna_bases(const std::string &dna_str);

    Pattern encode_dna_bases(const GenomicRead &read_sequence);

    /**
     * Integer encode (range: 1-4) a dna base character.
     */
    Base encode_dna_base(const char &base_str);

    std::string full_path(const std::string &gram_dirpath,
                          const std::string &file_name);
}

#endif //GRAMTOOLS_UTILS_HPP
