#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <iostream>


#include <boost/functional/hash.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "variants.hpp"
#include "ranks.hpp"

using namespace sdsl;
using namespace std;


#ifndef GRAMTOOLS_BWT_SEARCH_H
#define GRAMTOOLS_BWT_SEARCH_H


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

void calc_kmer_matches(const FM_Index &fm_index, int k,
                       sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>> &kmer_idx,
                       sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>> &kmer_idx_rev,
                       sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> &kmer_sites,
                       std::vector<int> &mask_a, uint64_t maxx, sequence_set<std::vector<uint8_t>> &kmers_in_ref,
                       std::vector<std::vector<uint8_t>> &kmerfile, const VariantMarkers &variants,
                       const DNA_Rank &rank_all);

void precalc_ranks(const FM_Index &fm_index, const DNA_Rank& rank_all);

std::pair<uint64_t, uint64_t> bidir_search(const uint8_t next_char,
                                           const std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                                           const std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                                           const FM_Index &fm_index,
                                           const DNA_Rank &rank_all);

std::pair<uint32_t, std::vector<int>> get_location(const FM_Index &fm_index,
                                                   const uint64_t marker_idx, const uint64_t marker,
                                                   const bool last, std::vector<int> &allele,
                                                   const std::vector<int> &mask_a);

bool skip(const FM_Index &fm_index,
          uint64_t& left, uint64_t& right,
          uint64_t& left_rev, uint64_t& right_rev,
          uint64_t num, uint64_t maxx);

void bidir_search_bwd(const FM_Index &fm_index,
                      uint64_t left, uint64_t right,
                      uint64_t left_rev, uint64_t right_rev, // not used in bwd
                      const std::vector<uint8_t>::iterator fasta_pattern_begin,
                      const std::vector<uint8_t>::iterator fasta_pattern_end,
                      std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                      std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev, // not used in bwd
                      std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                      const std::vector<int> &mask_a, const uint64_t maxx, bool &first_del,
                      const bool kmer_precalc_done, const VariantMarkers &variants,
                      const DNA_Rank &rank_all);

void bidir_search_bwd(const FM_Index &fm_index,
                      uint64_t left, uint64_t right,
                      uint64_t left_rev, uint64_t right_rev, // not used in bwd
                      const std::vector<uint8_t>::iterator fasta_pattern_begin,
                      const std::vector<uint8_t>::iterator fasta_pattern_end,
                      std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                      std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev, // not used in bwd
                      std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                      const std::vector<int> &mask_a,
                      const uint64_t maxx,
                      bool &first_del,
                      const bool kmer_precalc_done,
                      const VariantMarkers &variants,
                      const DNA_Rank &rank_all,
                      const int thread_id);

std::vector<uint8_t>::iterator bidir_search_fwd(const FM_Index &fm_index,
                                                uint64_t left, uint64_t right,
                                                uint64_t left_rev, uint64_t right_rev,
                                                std::vector<uint8_t>::iterator pat_begin,
                                                std::vector<uint8_t>::iterator pat_end,
                                                std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                                                std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                                std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                                                std::vector<int> &mask_a, uint64_t maxx, bool &first_del,
                                                bool kmer_precalc_done, const DNA_Rank &rank_all);

#endif //GRAMTOOLS_BWT_SEARCH_H
