#include "ranks.hpp"
#include "kmers.hpp"


#ifndef GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP
#define GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP


void bidir_search_bwd(SA_Intervals &sa_intervals, uint64_t left, uint64_t right, Sites &sites, bool &delete_first_interval,
                      const std::vector<uint8_t>::iterator fasta_pattern_begin,
                      const std::vector<uint8_t>::iterator fasta_pattern_end, const std::vector<int> &allele_mask,
                      const uint64_t maxx, const bool kmer_precalc_done, const DNA_Rank &rank_all, const FM_Index &fm_index,
                      const int thread_id=0);

void process_matches_overlapping_variants(SA_Intervals &sa_intervals, SA_Intervals::iterator &sa_interval_it, Sites &sites,
                                          Sites::iterator &sites_it, const bool delete_first, const uint64_t maxx,
                                          const std::vector<int> &allele_mask, const FM_Index &fm_index,
                                          const int thread_id);

bool match_next_charecter(bool delete_first_interval, SA_Intervals &sa_intervals, SA_Intervals::iterator &sa_interval_it,
                          Sites &sites, Sites::iterator &sites_it, const uint8_t next_char, const DNA_Rank &rank_all,
                          const FM_Index &fm_index, const int thread_id);

bool match_next_charecter(bool delete_first_interval, SA_Intervals &sa_intervals, SA_Intervals::iterator &sa_interval_it,
                          Sites &sites, Sites::iterator &sites_it, const uint8_t next_char, const DNA_Rank &rank_all,
                          const FM_Index &fm_index, const int thread_id);

void add_sa_interval_for_skip(uint64_t previous_marker, SA_Intervals::iterator &sa_interval_it, uint64_t &last_begin,
                              bool &second_to_last, const uint64_t marker_idx, const uint64_t marker,
                              uint64_t &left_new, uint64_t &right_new, bool &ignore);

void update_sites_crossed_by_reads(SA_Intervals &sa_intervals, SA_Intervals::iterator &sa_interval_it, uint64_t left_new,
                                   uint64_t right_new, Sites &sites, Sites::iterator &sites_it, bool &second_to_last,
                                   bool ignore, bool last, uint64_t &last_begin, const std::vector<int> &allele_mask,
                                   const bool &delete_first, const uint64_t marker, const uint64_t marker_idx,
                                   const FM_Index &fm_index);

std::pair<uint64_t, std::vector<int>> get_variant_site_edge(std::vector<int> &allele,
                                                            const uint64_t marker,
                                                            const uint64_t marker_idx,
                                                            const std::vector<int> &allele_mask,
                                                            const bool last,
                                                            const FM_Index &fm_index);


#endif //GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP
