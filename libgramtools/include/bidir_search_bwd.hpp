#include "ranks.hpp"
#include "kmers.hpp"
#include "prg.hpp"


#ifndef GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP
#define GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP

void bidir_search_bwd(SA_Intervals &sa_intervals,
                      Sites &sites,
                      bool &delete_first_interval,
                      const bool kmer_precalc_done,
                      const std::vector<uint8_t>::const_iterator read_begin,
                      const std::vector<uint8_t>::const_iterator read_end,
                      const PRG_Info &prg_info);

bool reduce_search_scope(const uint8_t read_char, SA_Intervals &sa_intervals, Sites &sites,
                         const bool delete_first_interval,
                         const bool kmer_precalc_done, const bool read_char_is_last, const PRG_Info &prg_info);

void process_reads_overlapping_variants(SA_Intervals &sa_intervals, SA_Interval &sa_interval, Sites &sites, Site &site,
                                        const bool delete_first, const PRG_Info &prg_info);

bool match_next_charecter(const uint8_t next_char, bool delete_first_interval, SA_Intervals &sa_intervals,
                          SA_Intervals::iterator &sa_intervals_it, Sites &sites, Sites::iterator &sites_it,
                          const DNA_Rank &rank_all, const FM_Index &fm_index);

bool match_next_charecter(const uint8_t next_char, bool delete_first_interval, SA_Intervals &sa_intervals,
                          SA_Intervals::iterator &sa_intervals_it, Sites &sites, Sites::iterator &sites_it,
                          const DNA_Rank &rank_all, const FM_Index &fm_index);

void add_sa_interval_for_skip(uint64_t &left_new, uint64_t &right_new, uint64_t &last_begin, bool &second_to_last,
                              bool &ignore, const uint64_t marker, const uint64_t marker_idx,
                              const uint64_t previous_marker);

bool update_sites_crossed_by_reads(SA_Intervals &sa_intervals, Sites &sites, const uint64_t left_new,
                                   const uint64_t right_new, const bool at_first_sa_interval, const bool ignore,
                                   const bool last, const uint64_t last_begin, const bool second_to_last,
                                   const bool delete_first, const uint64_t marker, const uint64_t marker_idx,
                                   const PRG_Info &prg_info);

VariantSite get_variant_site_edge(std::vector<int> &allele, const uint64_t marker, const uint64_t marker_idx, const bool last,
                                  const PRG_Info &prg_info);


#endif //GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP
