#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "prg.hpp"
#include "parameters.hpp"
#include "utils.hpp"


#ifndef GRAMTOOLS_KMERS_HPP
#define GRAMTOOLS_KMERS_HPP

template<typename SEQUENCE>
struct sequence_hash {
    std::size_t operator()(const SEQUENCE &seq) const {
        std::size_t hash = 0;
        boost::hash_range(hash, seq.begin(), seq.end());
        return hash;
    }
};

template<typename SEQUENCE>
struct sequence_ordering_condition {
    std::size_t operator()(const SEQUENCE &lhs, const SEQUENCE &rhs) const {
        int64_t i = 0;
        while (i >= 0) {
            if (lhs[i] < rhs[i]) {
                return true;
            } else if (lhs[i] == rhs[i]) {
                i++;
            } else if (lhs[i] > rhs[i]) {
                return false;
            }
            if (i == lhs.size())
                break;
        }
        return false;
    }
};

template<typename SEQUENCE>
using unordered_vector_set = std::unordered_set<SEQUENCE, sequence_hash<SEQUENCE>>;

template<typename SEQUENCE>
using ordered_vector_set = std::set<SEQUENCE, sequence_ordering_condition<SEQUENCE>>;

using PrgIndexRange = std::pair<uint64_t, uint64_t>;
using KmerSuffixDiffs = std::vector<sdsl::int_vector<8>>;

using Sequence = std::vector<Base>;
using SequencesList = std::vector<Sequence>;

std::vector<PrgIndexRange> get_boundary_marker_indexes(const PRG_Info &prg_info);

std::vector<PrgIndexRange> get_kmer_region_ranges(std::vector<PrgIndexRange> &boundary_marker_indexes,
                                                  const uint64_t max_read_size,
                                                  const PRG_Info &prg_info);

uint64_t find_site_end_boundary(const uint64_t within_site_index,
                                const PRG_Info &prg_info);

SequencesList get_site_ordered_alleles(const uint64_t within_site_index,
                                       const PRG_Info &prg_info);

std::list<uint64_t> sites_inrange_left(const uint64_t outside_site_start_index,
                                       const uint64_t kmer_size,
                                       const PRG_Info &prg_info);

std::pair<uint64_t, uint64_t> get_nonvariant_region(const uint64_t site_end_boundary_index,
                                                    const PRG_Info &prg_info);

Sequence right_intersite_nonvariant_region(const uint64_t site_end_boundary_index,
                                           const PRG_Info &prg_info);

unordered_vector_set<Sequence> get_reverse_kmers_from_region(const PrgIndexRange &kmer_region_range,
                                                             const uint64_t kmer_size,
                                                             const PRG_Info &prg_info);

uint64_t find_site_start_boundary(const uint64_t &end_boundary_index,
                                  const PRG_Info &prg_info);

std::list<SequencesList> get_kmer_size_region_parts(const uint64_t &current_range_end_index,
                                                    const std::list<uint64_t> &inrange_sites,
                                                    const uint64_t kmer_size,
                                                    const PRG_Info &prg_info);

SequencesList get_paths_from_parts(const std::list<SequencesList> &region_parts);

bool update_allele_index_path(std::vector<uint64_t> &allele_current_index,
                              const std::vector<uint64_t> &number_of_alleles);

unordered_vector_set<Sequence> get_reverse_kmers_from_path(const Sequence &path,
                                                           const uint64_t kmer_size);

unordered_vector_set<Sequence> extract_variant_reverse_kmers(uint64_t &current_range_end_index,
                                                             const std::list<uint64_t> &inrange_sites,
                                                             const uint64_t kmer_size,
                                                             const PRG_Info &prg_info);

std::vector<PrgIndexRange> combine_overlapping_regions(const std::vector<PrgIndexRange> &kmer_region_ranges);

std::vector<Sequence> reverse_kmers_inplace(const ordered_vector_set<Sequence> &reverse_kmers);

std::vector<Sequence> get_prefix_diffs(const std::vector<Sequence> &kmers);

ordered_vector_set<Sequence> get_all_reverse_kmers(const Parameters &parameters,
                                                   const PRG_Info &prg_info);

std::vector<Sequence> get_kmer_prefix_diffs(const Parameters &parameters,
                                            const PRG_Info &prg_info);

std::vector<Sequence> get_all_ordered_kmers(const Parameters &parameters,
                                            const PRG_Info &prg_info);

#endif //GRAMTOOLS_KMERS_HPP
