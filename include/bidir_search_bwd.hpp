#include "ranks.hpp"


#ifndef GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP
#define GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP


void bidir_search_bwd(std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                      std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev,
                      uint64_t left, uint64_t right,
                      uint64_t left_rev, uint64_t right_rev,
                      std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                      bool &delete_first_interval,
                      const std::vector<uint8_t>::iterator fasta_pattern_begin,
                      const std::vector<uint8_t>::iterator fasta_pattern_end,
                      const std::vector<int> &mask_a,
                      const uint64_t maxx,
                      const bool kmer_precalc_done,
                      const DNA_Rank &rank_all,
                      const FM_Index &fm_index,
                      const int thread_id);

void update_sites_crossed_by_reads(std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                                   std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                   std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                                   std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                                   uint64_t left_new, uint64_t right_new, uint64_t left_rev_new, uint64_t right_rev_new,
                                   std::vector<int> &allele_empty,
                                   std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                                   std::list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                                   bool &second_to_last, bool ignore, bool last, uint64_t &last_begin,
                                   const std::vector<int> &mask_a, const bool &first_del, const uint64_t marker,
                                   const uint64_t marker_idx, const FM_Index &fm_index);

bool match_next_charecter(bool delete_first_interval,
                          std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                          std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev,
                          std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                          std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                          std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                          std::list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                          const uint8_t next_char,
                          const DNA_Rank &rank_all,
                          const FM_Index &fm_index,
                          const int thread_id);

void add_sa_interval_for_skip(uint64_t previous_marker,
                              std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                              std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                              uint64_t &last_begin, bool &second_to_last,
                              const uint64_t marker_idx, const uint64_t marker,
                              uint64_t &left_new, uint64_t &right_new,
                              uint64_t &left_rev_new, uint64_t &right_rev_new,
                              bool &ignore);

void process_matches_overlapping_variants(std::vector<int> &allele_empty,
                                          std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                                          std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                                          std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                                          std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                          std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                                          std::list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                                          const bool first_del,
                                          const uint64_t maxx,
                                          const std::vector<int> &mask_a,
                                          const FM_Index &fm_index,
                                          const int thread_id);

std::pair<uint64_t, std::vector<int>> get_variant_site_edge(std::vector<int> &allele,
                                                            const uint64_t marker,
                                                            const uint64_t marker_idx,
                                                            const std::vector<int> &mask_a,
                                                            const bool last,
                                                            const FM_Index &fm_index);


#endif //GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP
