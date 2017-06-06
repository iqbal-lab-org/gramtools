#ifndef __BWT_SEARCH_H_INCLUDED__
#define __BWT_SEARCH_H_INCLUDED__

#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <boost/functional/hash.hpp>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <iostream>

using namespace sdsl;
using namespace std;


using CSA = csa_wt<wt_int<bit_vector, rank_support_v5<>>, 2, 16777216>;

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

CSA csa_constr(std::string fname, std::string int_al_fname,
               std::string memory_log_fname, std::string csa_file,
               bool fwd, bool verbose);

void precalc_kmer_matches(CSA &csa, int k,
                          sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>> &kmer_idx,
                          sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>> &kmer_idx_rev,
                          sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> &kmer_sites,
                          std::vector<int> &mask_a, uint64_t maxx, sequence_set<std::vector<uint8_t>> &kmers_in_ref,
                          std::vector<std::vector<uint8_t>> &kmerfile);

uint64_t bidir_search(CSA &csa,
                      uint64_t &left, uint64_t &right,
                      uint64_t &left_rev, uint64_t &right_rev,
                      uint8_t c);


std::pair<uint32_t, std::vector<int>> get_location(CSA &csa,
                                                   uint64_t num_idx,
                                                   uint32_t num, bool last,
                                                   std::vector<int> &allele,
                                                   std::vector<int> &mask_a);

bool skip(CSA &csa,
          uint64_t &left, uint64_t &right,
          uint64_t &left_rev, uint64_t &right_rev,
          uint32_t num, uint64_t maxx);

std::vector<uint8_t>::iterator bidir_search_bwd(CSA &csa,
                                                uint64_t left, uint64_t right,
                                                uint64_t left_rev, uint64_t right_rev,
                                                std::vector<uint8_t>::iterator pat_begin,
                                                std::vector<uint8_t>::iterator pat_end,
                                                std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                                                std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                                std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                                                std::vector<int> &mask_a, uint64_t maxx, bool &first_del,
                                                bool kmer_precalc_done);


std::vector<uint8_t>::iterator bidir_search_fwd(CSA &csa,
                                                uint64_t left, uint64_t right,
                                                uint64_t left_rev, uint64_t right_rev,
                                                std::vector<uint8_t>::iterator pat_begin,
                                                std::vector<uint8_t>::iterator pat_end,
                                                std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                                                std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                                std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                                                std::vector<int> &mask_a, uint64_t maxx, bool &first_del,
                                                bool kmer_precalc_done);

#endif
