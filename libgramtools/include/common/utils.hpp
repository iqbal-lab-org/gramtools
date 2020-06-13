/** @file
 * Defines a set of data structures and functions common to the backend.
 */
#ifndef GRAMTOOLS_UTILS_HPP
#define GRAMTOOLS_UTILS_HPP

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>
#include "sequence_read/seqread.hpp"

namespace gram {
using int_Base = uint8_t; /**< nucleotide represented as byte-sized integer */
using Sequence =
    std::vector<int_Base>; /** A string of nucleotides is represented as a
                              vector of `Base`s. */
using Sequences = std::vector<Sequence>;

/******************************************
 * Characters to integers and vice-versa **
 ******************************************/
struct EncodeResult {
  bool is_dna;
  uint32_t character;
};

/**
 * Encode a single character read from the linearised prg as an integer.
 * Use `EncodeResult` object to additionally store if the encoded character is
 * DNA or not.
 * @see EncodeResult()
 */
EncodeResult encode_char(const char &c);

/**
 * Encode dna character as integer (range: 1-4)
 */
int_Base encode_dna_base(const char &base_str);

/**
 * Decode integer into dna base, as std::string.
 */
std::string decode_dna_base(const int_Base &base);

Sequence encode_dna_bases(const std::string &dna_str);

Sequence encode_dna_bases(const GenomicRead &read_sequence);

/************
 * Hashing **
 ************/
template <typename SEQUENCE>
struct sequence_hash {
  std::size_t operator()(const SEQUENCE &seq) const {
    std::size_t hash = 0;  // Used as an initial seed
    boost::hash_range(hash, seq.begin(), seq.end());
    return hash;
  }
};

template <typename SEQUENCE, typename T>
using SequenceHashMap =
    std::unordered_map<SEQUENCE, T, sequence_hash<SEQUENCE>>;

template <typename PAIR, typename T>
using PairHashMap = std::unordered_map<PAIR, T, boost::hash<PAIR>>;

template <typename T>
using HashSet = std::unordered_set<T, boost::hash<T>>;

}  // namespace gram

#endif  // GRAMTOOLS_UTILS_HPP
