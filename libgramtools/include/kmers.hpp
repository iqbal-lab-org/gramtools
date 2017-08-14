#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include "fm_index.hpp"
#include "ranks.hpp"


#ifndef GRAMTOOLS_KMERS_HPP
#define GRAMTOOLS_KMERS_HPP


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

using Allele = std::vector<int>;
using VariantSiteMarker = uint64_t;
using VariantSite = std::pair<VariantSiteMarker, Allele>;
using Site = std::vector<VariantSite>;
using Sites = std::list<Site>;

using SA_Interval = std::pair<uint64_t, uint64_t>;
using SA_Intervals = std::list<SA_Interval>;

using Kmer = std::vector<uint8_t>;
using Kmers = std::vector<Kmer>;
using KmerIdx = sequence_map<Kmer, SA_Intervals>;
using KmerSites = sequence_map<Kmer, Sites>;
using KmersRef = sequence_set<Kmer>;

inline bool fexists(const std::string &name);
static inline std::string &ltrim(std::string &s);
static inline std::string &rtrim(std::string &s);
static inline std::string &trim(std::string &s);
std::vector<std::string> split(const std::string &cad, const std::string &delim);

struct ThreadData {
    KmerIdx *kmer_idx;
    KmerSites *kmer_sites;
    KmersRef *kmers_in_ref;
    std::vector<std::vector<uint8_t>> *kmers;
    int thread_id;
    const FM_Index *fm_index;
    const DNA_Rank *rank_all;
    const std::vector<int> *allele_mask;
    uint64_t maxx;
};

void *worker(void *st);

struct KmersData {
    KmerIdx index;
    KmerSites sites;
    KmersRef in_reference;
};

KmersData read_encoded_kmers(const std::string &fil);

KmersData get_kmers(const std::string &kmer_fname, const std::vector<int> &allele_mask, const uint64_t maxx,
                    const DNA_Rank &rank_all, const FM_Index &fm_index);

void gen_precalc_kmers(const FM_Index &fm_index,
                       const std::vector<int> &allele_mask,
                       const std::string &kmer_fname,
                       const uint64_t maxx,
                       const int k,
                       const DNA_Rank &rank_all);

void calc_kmer_matches(KmerIdx &kmer_idx, KmerSites &kmer_sites, KmersRef &kmers_in_ref, Kmers &kmers, const uint64_t maxx,
                       const std::vector<int> &allele_mask, const DNA_Rank &rank_all, const FM_Index &fm_index,
                       const int thread_id);


#endif //GRAMTOOLS_KMERS_HPP
