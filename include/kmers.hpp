#include "bwt_search.h"
#include "variants.hpp"


#ifndef GRAMTOOLS_KMERS_HPP
#define GRAMTOOLS_KMERS_HPP


using KmerIdx = sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>>;
using KmerSites = sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>;
using KmersRef = sequence_set<std::vector<uint8_t>>;

inline bool fexists(const std::string &name);

//trim from start
static inline std::string &ltrim(std::string &s);

// trim from end
static inline std::string &rtrim(std::string &s);

// trim from both ends
static inline std::string &trim(std::string &s);

std::vector<std::string> split(std::string cad, std::string delim);

struct thread_data {
    VariantMarkers *variants;
    csa_wt<wt_int<bit_vector, rank_support_v5<>>, 2, 16777216> *csa;
    int k;
    KmerIdx *kmer_idx, *kmer_idx_rev;
    KmerSites *kmer_sites;
    std::vector<int> *mask_a;
    uint64_t maxx;
    KmersRef *kmers_in_ref;
    std::vector<std::vector<uint8_t>> *kmers;
    int thread_id;
};

void *worker(void *st);

void gen_precalc_kmers(csa_wt<wt_int<bit_vector, rank_support_v5<>>, 2, 16777216> &csa,
                       std::vector<int> &mask_a, std::string kmer_fname,
                       uint64_t maxx, int k, VariantMarkers &variants);

void read_precalc_kmers(std::string fil, KmerIdx &kmer_idx,
                        KmerIdx &kmer_idx_rev, KmerSites &kmer_sites,
                        KmersRef &kmers_in_ref);

struct KmersData {
    KmerIdx index, index_reverse;
    KmerSites sites;
    KmersRef in_reference;
};

KmersData get_kmers(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> &csa,
                    std::vector<int> &mask_a, std::string kmer_fname,
                    uint64_t maxx, int k, VariantMarkers &variants);


#endif //GRAMTOOLS_KMERS_HPP
