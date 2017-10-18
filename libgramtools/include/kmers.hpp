#include "utils.hpp"
#include "prg.hpp"
#include "fm_index.hpp"
#include "ranks.hpp"
#include "search.hpp"
#include "kmer_index.hpp"


#ifndef GRAMTOOLS_KMERS_HPP
#define GRAMTOOLS_KMERS_HPP

struct CacheElement {
    SearchStates search_states;
    Base base = 0;
};

using KmerIndexCache = std::list<CacheElement>;

std::ostream &operator<< (std::ostream &os, const KmerIndexCache &cache);

Base encode_dna_base(const char &base_str);

std::vector<Base> encode_dna_bases(const std::string &dna_str);

std::string dump_kmer(const Pattern &kmer);

std::string dump_sa_intervals(const SA_Intervals &sa_intervals);

std::string dump_crosses_marker_flag(const Pattern &kmer,
                                     const NonSiteCrossingKmers &nonvar_kmers);

std::string dump_variant_site_paths(const Pattern &kmer, const KmerVariantSitePaths &kmer_sites);

std::string dump_kmer_index_entry(const Pattern &kmer, const SA_Intervals &sa_intervals, const KmerVariantSitePaths &kmer_sites,
                                  const NonSiteCrossingKmers &nonvar_kmers);

void dump_kmer_index(std::ofstream &precalc_file, const KmerIndex &kmer_index);

KmerIndex index_kmers(const Patterns &kmers, const int kmer_size, const PRG_Info &prg_info);

inline bool file_exists(const std::string &name);

static inline std::string &left_trim(std::string &s);

static inline std::string &right_trim(std::string &s);

static inline std::string &trim(std::string &s);

std::vector<std::string> split(const std::string &cad, const std::string &delim);

bool parse_crosses_marker_flag(const std::string &in_reference_flag_str);

Pattern parse_encoded_kmer(const std::string &encoded_kmer_str);

SA_Intervals parse_sa_intervals(const std::string &full_sa_intervals_str);

VariantSitePath parse_variant_site_path(const std::string &sites_part_str);

void parse_kmer_index_entry(KmerIndex &kmers, const std::string &line);

KmerIndex load_kmer_index(const std::string &encoded_kmers_fname);

KmerIndex get_kmer_index(const std::string &kmer_fname, const int kmer_size, const PRG_Info &prg_info);

#endif //GRAMTOOLS_KMERS_HPP
