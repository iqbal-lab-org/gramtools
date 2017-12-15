#include "utils.hpp"
#include "parameters.hpp"
#include "prg.hpp"
#include "fm_index.hpp"
#include "ranks.hpp"
#include "search_types.hpp"
#include "kmer_index_types.hpp"
#include "kmers.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_HPP
#define GRAMTOOLS_KMER_INDEX_HPP

std::ostream &operator<<(std::ostream &os, const KmerIndexCache &cache);

Base encode_dna_base(const char &base_str);

std::string dump_kmer(const Pattern &kmer);

std::string dump_sa_intervals(const SearchStates &search_states);

std::string dump_variant_site_paths(const SearchStates &kmer_sites);

std::string dump_kmer_index_entry(const Pattern &kmer, const SearchStates &search_states);

void dump_kmer_index(std::ofstream &kmer_index_file, const KmerIndex &kmer_index);

KmerIndex index_kmers(const Patterns &kmers, const int kmer_size, const PRG_Info &prg_info);

inline bool file_exists(const std::string &name);

static inline std::string &left_trim(std::string &s);

static inline std::string &right_trim(std::string &s);

static inline std::string &trim(std::string &s);

std::vector<std::string> split(const std::string &cad, const std::string &delim);

bool parse_crosses_marker_flag(const std::string &in_reference_flag_str);

Pattern parse_encoded_kmer(const std::string &encoded_kmer_str);

std::vector<SA_Interval> parse_sa_intervals(const std::string &full_sa_intervals_str);

VariantSitePath parse_variant_site_path(const std::string &sites_part_str);

void parse_kmer_index_entry(KmerIndex &kmers, const std::string &line);

KmerIndex load_kmer_index(const Parameters &parameters);

void generate_kmer_index(const Parameters &parameters,
                         const PRG_Info &prg_info);

KmerIndex get_kmer_index(const std::string &kmer_fname, const int kmer_size, const PRG_Info &prg_info);

#endif //GRAMTOOLS_KMER_INDEX_HPP
