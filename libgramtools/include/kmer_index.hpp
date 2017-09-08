#include <sdsl/suffix_trees.hpp>
#include "fm_index.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_HPP
#define GRAMTOOLS_KMER_INDEX_HPP

// using WavletTree = sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>>;
// using FM_Index = sdsl::csa_wt<WavletTree, 2, 16777216>;
//using KmerIndex = sdsl::cst_sada<>;

void generate_kmer_index(const Kmers &kmers, const uint64_t max_alphabet_number, const std::vector<int> &allele_mask,
                         const DNA_Rank &rank_all, const FM_Index &fm_index);


#endif //GRAMTOOLS_KMER_INDEX_HPP
