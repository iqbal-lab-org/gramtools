#include "bwt_search.hpp"


#ifndef GRAMTOOLS_KMERS_HPP
#define GRAMTOOLS_KMERS_HPP


inline bool fexists(const std::string &name);

//trim from start
static inline std::string &ltrim(std::string &s);

// trim from end
static inline std::string &rtrim(std::string &s);

// trim from both ends
static inline std::string &trim(std::string &s);

std::vector<std::string> split(const std::string &cad, const std::string &delim);

using KmerIdx = sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>>;
using KmerSites = sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint64_t, std::vector<int>>>>>;
using KmersRef = sequence_set<std::vector<uint8_t>>;

struct ThreadData {
    KmerIdx *kmer_idx, *kmer_idx_rev;
    KmerSites *kmer_sites;
    KmersRef *kmers_in_ref;
    std::vector<std::vector<uint8_t>> *kmers;
    int thread_id;
    const FM_Index *fm_index;
    const DNA_Rank *rank_all;
    int k;
    const std::vector<int> *mask_a;
    uint64_t maxx;
};

void *worker(void *st);

struct KmersData {
    KmerIdx index, index_reverse;
    KmerSites sites;
    KmersRef in_reference;
};

KmersData read_encoded_kmers(const std::string &fil);

KmersData get_kmers(const std::string &kmer_fname, const int k, const std::vector<int> &mask_a, const uint64_t maxx,
                    const DNA_Rank &rank_all, const FM_Index &fm_index);

void gen_precalc_kmers(const FM_Index &fm_index,
                       const std::vector<int> &mask_a,
                       const std::string &kmer_fname,
                       const uint64_t maxx,
                       const int k,
                       const DNA_Rank &rank_all);

void calc_kmer_matches(KmerIdx &kmer_idx,
                       KmerIdx &kmer_idx_rev,
                       KmerSites &kmer_sites,
                       sequence_set<std::vector<uint8_t>> &kmers_in_ref,
                       std::vector<std::vector<uint8_t>> &kmers,
                       const FM_Index &fm_index,
                       const DNA_Rank &rank_all,
                       const std::vector<int> &mask_a,
                       const int k,
                       const uint64_t maxx);


#endif //GRAMTOOLS_KMERS_HPP
