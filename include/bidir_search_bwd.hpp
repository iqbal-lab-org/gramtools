#ifndef GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP
#define GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP


void update_sites_crossed_by_reads(const FM_Index &fm_index, list<pair<uint64_t, uint64_t>> &sa_intervals,
                                   list<pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                   list<vector<pair<uint32_t, vector<int>>>> &sites, const vector<int> &mask_a,
                                   const bool &first_del, vector<int> &allele_empty, bool &second_to_last,
                                   uint64_t marker_idx, uint64_t marker, uint64_t left_new, uint64_t right_new,
                                   uint64_t left_rev_new, uint64_t right_rev_new, bool ignore, bool last,
                                   list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                                   list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                                   list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                                   uint64_t &last_begin, MarkerPositions &marker_it);

bool match_next_charecter(const FM_Index &fm_index, list<pair<uint64_t, uint64_t>> &sa_intervals,
                          list<pair<uint64_t, uint64_t>> &sa_intervals_rev,
                          list<vector<pair<uint32_t, vector<int>>>> &sites,
                          uint8_t fasta_char,
                          list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                          list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                          list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                          bool first_del);

void add_sa_interval_for_skip(uint64_t previous_marker,
                              list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                              list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                              uint64_t &last_begin, bool &second_to_last,
                              MarkerPositions &marker_it,
                              uint64_t &marker_idx, uint64_t &marker, uint64_t &left_new, uint64_t &right_new,
                              uint64_t &left_rev_new, uint64_t &right_rev_new, bool &ignore);

void process_matches_overlapping_variants(const FM_Index &fm_index, const vector<int> &mask_a, const uint64_t maxx,
                                          const bool &first_del, const VariantMarkers &variants,
                                          vector<int> &allele_empty,
                                          list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                                          list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                                          list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                                          list<pair<uint64_t, uint64_t>> &sa_intervals,
                                          list<pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                          list<vector<pair<uint32_t, vector<int>>>> &sites);


#endif //GRAMTOOLS_BIDIR_SEARCH_BWD_HPP_HPP
