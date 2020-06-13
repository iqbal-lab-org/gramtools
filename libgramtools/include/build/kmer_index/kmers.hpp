/**
 * @file
  * Routines for producing all required kmers, in sorted order.
 /******************
 /* * Generalities *
 /******************
 * Kmers are produces either from:
 * * The prg, extracting only kmers overlapping variant sites
 * * All possible kmers
 *
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
 *
 /******************************
 /* * Kmer extraction from prg *
 /******************************
 * The general approach is as follows:
 * * Find all variant site start-end positions in the prg.
 * * Extend them to the right up to the maximum read size, such that all kmers
 in a read which could end in a variant site are included.
 * * Combine overlapping above regions to avoid redundancy.
 * * For each position of each region:
 *      * Find all variant sites to the left within range
 *      * Enumerate all possible paths in the prg going through these sites
 *      * Extract all kmers of the given (user-defined) size in each path and
 add them to the set of kmers to index.
 */
#include <unordered_map>
#include <unordered_set>

#include <boost/functional/hash.hpp>

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

using PrgIndexRange = std::pair<uint64_t, uint64_t>;
using KmerSuffixDiffs = std::vector<sdsl::int_vector<8>>;

/**
 * Computes the start and end indices of all variant site markers in a given
 * prg.
 */
std::vector<PrgIndexRange> get_boundary_marker_indexes(
    const PRG_Info &prg_info);

/**
 * Extend each variant site region in the prg to include read-reachable regions.
 * The rationale is to index only kmers whose extension backward in the prg by
 * mapping
 * **can** overlap a variant site.
 * @param max_read_size the maximum size of a read to map to the prg.
 */
std::vector<PrgIndexRange> get_kmer_region_ranges(
    std::vector<PrgIndexRange> &boundary_marker_indexes,
    const uint64_t &max_read_size, const PRG_Info &prg_info);

/**
 * Finds the index, in the prg, of the end boundary of a variant site.
 * @note Passing in an end index returns the end index itself.
 */
uint64_t find_site_end_boundary(const uint64_t &within_site_index,
                                const PRG_Info &prg_info);

/**
 * Extract all alleles of a variant site.
 * @param within_site_index
 * @return A vector of alleles. Alleles are ordered by appearance in the prg.
 */
Sequences get_site_ordered_alleles(const uint64_t &within_site_index,
                                   const PRG_Info &prg_info);

/**
 * Record which variant sites are potentially reachable by a kmer.
 * From the starting index, sequence is consumed up to the kmer size.
 * @param outside_site_start_index the starting position. This can only be
 * outside of a variant site.
 * @return list of the end indices in the prg of the variant sites in range of
 * the kmer.
 */
std::list<uint64_t> sites_inrange_left(const uint64_t outside_site_start_index,
                                       const uint64_t kmer_size,
                                       const PRG_Info &prg_info);

std::pair<uint64_t, uint64_t> get_nonvariant_region(
    const uint64_t &site_end_boundary_index, const PRG_Info &prg_info);

Sequence right_intersite_nonvariant_region(
    const uint64_t &site_end_boundary_index, const PRG_Info &prg_info);

/** Build a set of kmers to index from each index inside a `kmer_region_range`.
 * @param kmer_region_range Set of contiguous indices into the prg from which
 * sets of kmers to be indexed are extracted ony by one.
 */
unordered_vector_set<Sequence> get_region_range_reverse_kmers(
    const PrgIndexRange &kmer_region_range, const uint64_t &kmer_size,
    const PRG_Info &prg_info);

/**
 * Find the start index of a variant site marker from its end index.
 */
uint64_t find_site_start_boundary(const uint64_t &end_boundary_index,
                                  const PRG_Info &prg_info);

/**
 * From a list of reachable variant sites, extract a set of parts (alleles and
 * non-variant regions). These parts are later combined to generate all the
 * kmers to index starting from a given prg index.
 */
std::list<Sequences> get_kmer_size_region_parts(
    const uint64_t &current_range_end_index,
    const std::list<uint64_t> &inrange_sites, const uint64_t kmer_size,
    const PRG_Info &prg_info);

/**
 * Extracts all kmers to index from the `region_parts` that can be traversed.
 * @param region_parts The variant and non-variant regions within reach of a
 * starting position in the prg.
 * @return all unique kmers extracted from all enumerated path through the
 * `region_parts`.
 */
unordered_vector_set<Sequence> get_region_parts_reverse_kmers(
    const std::list<Sequences> &region_parts, const uint64_t &kmer_size);

/**
 * Increments a single allele index among all region parts.
 * This allows exhaustive run-through of all possible paths through variant
 * sites.
 */
bool update_allele_index_path(std::vector<uint64_t> &current_allele_index_path,
                              const std::vector<uint64_t> &parts_allele_counts);

/**
 * From a `path`, extract all kmers to index.
 * @param path A single path through the prg.
 * @return A (unique) set of kmers to index. The kmers are in reverse (right to
 * left) order.
 */
unordered_vector_set<Sequence> get_path_reverse_kmers(
    const Sequence &path, const uint64_t &kmer_size);

/**
 * Gets all unique kmers to index from a starting position in the prg.
 * The kmers are extracted from all possible paths crossing the variant sites in
 * range.
 * @param current_range_end_index the starting position in the prg.
 * @param inrange_sites the variant sites in range of the starting position.
 * @note the function updates `current_range_end_index` so that it is past the
 * leftmost site in `inrange_sites`.
 */
unordered_vector_set<Sequence> get_sites_reverse_kmers(
    uint64_t &current_range_end_index, const std::list<uint64_t> &inrange_sites,
    const uint64_t kmer_size, const PRG_Info &prg_info);

/**
 * Sort a set of kmer ranges (`gram::PrgIndexRange`s), and merge together those
 * that overlap. Produces maximally large, non-overlapping kmer ranges.
 * @note an overlap between two kmer regions means at least two distinct sites
 * can be mapped into by a single read.
 */
std::vector<PrgIndexRange> combine_overlapping_regions(
    const std::vector<PrgIndexRange> &kmer_region_ranges);

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
 * Extract kmers to index from a prg. Only kmers in the prg whose mapping can
 * overlap a variant site will get indexed.
 * @return all the kmers to index, in reverse sorted order. The kmers are
 * maintained in reverse (first kmer position == last position in prg) so that
 * they get inserted in sorted (dictionary) order.
 */
ordered_vector_set<Sequence> get_prg_reverse_kmers(
    BuildParams const &parameters, const PRG_Info &prg_info);

/**
 * High-level routine for extracting all kmers of interest and computing the
 * prefix differences.
 * @see gram::kmer_index::build()
 */
std::vector<Sequence> get_all_kmer_and_compute_prefix_diffs(
    BuildParams const &parameters, const PRG_Info &prg_info);

/**
 * Core routine for producing kmers to index.
 * If `all_kmers_flag` is unset (the default), only the relevant kmers are
 * produced; these are the kmers overlapping variant sites in the prg. If it is
 * set, produces all the kmers of given size.
 */
std::vector<Sequence> get_all_kmers(BuildParams const &parameters,
                                    const PRG_Info &prg_info);

/**
 * Generate all kmers of a given size, in order.
 * Order is a dictionary order ('1111'<'1121' < '1211' etc..)
 * Kmers themselves are produced in order ('1111' then '1112' etc..)
 * @see next_kmer()
 */
ordered_vector_set<Sequence> generate_all_kmers(const uint64_t &kmer_size);

}  // namespace gram

#endif  // GRAMTOOLS_KMERS_HPP
