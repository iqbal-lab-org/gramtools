#include "utils.hpp"
#include "kmer_index_types.hpp"
#include "search_types.hpp"


#ifndef GRAMTOOLS_SEARCH_HPP
#define GRAMTOOLS_SEARCH_HPP

#define USE_SKIP_OPTIMIZATION false


SearchStates search_read_backwards(const Pattern &read,
                                   const Pattern &kmer,
                                   const KmerIndex &kmer_index,
                                   const PRG_Info &prg_info);

SearchStates search_base_backwards(const Base &pattern_char,
                                   const SearchStates &search_states,
                                   const PRG_Info &prg_info);

using SaIndexRightOfMarker = uint64_t;
using MarkersSearchResult = std::pair<SaIndexRightOfMarker, Marker>;
using MarkersSearchResults = std::vector<MarkersSearchResult>;

MarkersSearchResults left_markers_search(const SearchState &search_state,
                                         const PRG_Info &prg_info);

SA_Interval base_next_sa_interval(const Marker current_char,
                                  const SA_Index &current_char_first_sa_index,
                                  const SA_Interval &current_sa_interval,
                                  const PRG_Info &prg_info);

void process_search_state_path_cache(SearchState &search_state);

SearchStates process_read_char_search_states(const Base &pattern_char,
                                             const SearchStates &old_search_states,
                                             const PRG_Info &prg_info);

SearchState search_skipping_marker(const SearchState &search_state,
                                   const Base &pattern_char,
                                   const PRG_Info &prg_info);

void set_state_skip_marker(SearchStates &search_states,
                           const PRG_Info &prg_info);

AlleleId get_allele_id(const SA_Index &allele_marker_sa_index,
                       const PRG_Info &prg_info);

SA_Interval get_allele_marker_sa_interval(const Marker &site_marker_char,
                                          const PRG_Info &prg_info);

SearchStates process_markers_search_states(const SearchStates &search_states,
                                           const PRG_Info &prg_info);

SearchStates process_markers_search_state(const SearchState &search_state,
                                          const PRG_Info &prg_info);

std::string serialize_search_state(const SearchState &search_state);

std::ostream &operator<< (std::ostream &os, const SearchState &search_state);

#endif //GRAMTOOLS_SEARCH_HPP
