/**
 * @file
 * Routines for producing all kmers to index, in sorted order.
 *
 * ** Output **
 * All possible kmers of a given size are produced.
 * It would be nice to only output only kmers overlapping variant sites
 * (and those at +- maximum read size).
 * Such a kmer list has exponential size in worst case.
 *
 * ** Implementation **
 * The kmers are sorted such that each kmer shares the largest possible suffix
 with its predecessor in the set.
 * This reduces the number of prg quasimappings to compute to a minimum. Indeed,
 all the search states associated with a
 * given kmer suffix are kept in cache and then extended to produce larger kmer
 suffixes.
 *
 * For example, kmers '1111' and '2111' (corresponding to 'aaaa' and 'caaa'
 respectively) are stored consecutively.
 * Then the set of `SearchStates` corresponding to the string '111' in the PRG
 is computed only once, and used to extend
 * to both '1111' and '2111' during variant aware backward searching.
 */
#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <unordered_set>

#include "build/parameters.hpp"
#include "common/utils.hpp"
#include "prg/prg_info.hpp"

#ifndef GRAMTOOLS_KMERS_HPP
#define GRAMTOOLS_KMERS_HPP

namespace gram {

template <typename SEQUENCE>
struct sequence_ordering_condition {
  std::size_t operator()(const SEQUENCE &lhs, const SEQUENCE &rhs) const {
    int64_t i = 0;
    while (i >= 0) {
      if (lhs[i] < rhs[i]) {
        return (std::size_t) true;
      } else if (lhs[i] == rhs[i]) {
        i++;
      } else if (lhs[i] > rhs[i]) {
        return (std::size_t) false;
      }
      if (i == lhs.size()) break;
    }
    return (std::size_t) false;
  }
};

template <typename SEQUENCE>
using unordered_vector_set =
    std::unordered_set<SEQUENCE, sequence_hash<SEQUENCE>>;

template <typename SEQUENCE>
using ordered_vector_set =
    std::set<SEQUENCE, sequence_ordering_condition<SEQUENCE>>;

/**
 * Converts an ordered_set of kmers into a vector of kmers in the reverse order.
 * The use of this is that the ordered_set when applied on kmers stored
 * right-to-left in the prg, naturally maximises the shared suffix between
 * consecutive entries. Reversing them so that they are left-to-right in the
 * prg, readies the kmers for the cached indexing process.
 * @see gram::get_prefix_diffs()
 */
std::vector<Sequence> reverse(
    const ordered_vector_set<Sequence> &reverse_kmers);

/**
 * Computes the minimal changes between kmers and their immediate predecessor in
 * the ordered set. The changes are listed left to right (prefix order).
 */
std::vector<Sequence> get_prefix_diffs(const std::vector<Sequence> &kmers);

/**
 * Generate all kmers of a given size, in dictionary order
 * ('1111'<'1121' < '1211' etc..).
 * @see next_kmer()
 */
ordered_vector_set<Sequence> generate_all_kmers(const uint64_t &kmers_size);

/**
 * Glue function
 * @see gram::get_all_kmer_and_compute_prefix_diffs()
 */
std::vector<Sequence> get_all_kmers(const uint64_t &kmers_size);

/**
 * Produces a list of all kmers to index and computes their prefix differences,
 * i.e. how different two consecutive kmers in the list are.
 * @see gram::kmer_index::build()
 */
std::vector<Sequence> get_all_kmer_and_compute_prefix_diffs(
    uint64_t const &kmers_size);

}  // namespace gram

#endif  // GRAMTOOLS_KMERS_HPP
