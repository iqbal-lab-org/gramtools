#include "sequence_read/seqread.hpp"

#include "bwt_search.hpp"
#include "parameters.hpp"
#include "masks.hpp"
#include "kmers.hpp"
#include "fm_index.hpp"
#include "ranks.hpp"


#ifndef GRAMTOOLS_MAP_HPP
#define GRAMTOOLS_MAP_HPP

void print_sa_interval(const SA_Intervals &sa_intervals);

void print_sites(const Sites &sites);

int quasimap_reads(KmerIndex &kmers, MasksParser &masks, const Parameters &params, const PRG_Info &prg_info);

std::vector<uint8_t> int_encode_read(const GenomicRead &read_sequence);

bool process_read(const std::vector<uint8_t> &encoded_read,
                  int &count_char_in_variant_site,
                  std::unordered_set<uint64_t> &repeats_variant_site_edge_markers,
                  KmerIndex &kmers,
                  MasksParser &masks,
                  const Parameters &params,
                  const DNA_Rank &rank_all,
                  const FM_Index &fm_index);

bool discard_read_check(const std::vector<uint8_t> &read_kmer_part, const KmerIndex &kmers);

bool quasimap_read(const std::vector<uint8_t> &read_kmer_part,
                   const std::vector<uint8_t> &encoded_read,
                   int &count_char_in_variant_site,
                   std::unordered_set<uint64_t> &repeats_variant_site_edge_markers,
                   KmerIndex &kmers,
                   MasksParser &masks,
                   const int kmer_size,
                   const DNA_Rank &rank_all,
                   const FM_Index &fm_index);

void populate_repeats_variant_edges(std::unordered_set<uint64_t> &repeats_variant_site_edge_markers,
                                    int &count_char_in_variant_site,
                                    const MasksParser &masks,
                                    const Sites &sites,
                                    const SA_Intervals &sa_intervals,
                                    const FM_Index &fm_index);

void update_site_sa_interval_coverage(MasksParser &masks, const Site &site, const bool first_site_empty,
                                      const SA_Interval &sa_interval, const bool is_first_sa_interval,
                                      const bool delete_first, const uint64_t total_num_sa_intervals,
                                      const int in_sites, const std::vector<uint8_t> &readin_integer_seq,
                                      const uint64_t repeats, const DNA_Rank &rank_all, const FM_Index &fm_index);

void output_allele_coverage(Parameters &params, MasksParser &masks);


#endif //GRAMTOOLS_MAP_HPP
