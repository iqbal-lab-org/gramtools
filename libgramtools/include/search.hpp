#include "utils.hpp"
#include "kmer_index_types.hpp"
#include "search_types.hpp"


#ifndef GRAMTOOLS_SEARCH_HPP
#define GRAMTOOLS_SEARCH_HPP

std::string serialize_search_state(const SearchState &search_state);

std::ostream &operator<< (std::ostream &os, const SearchState &search_state);

SearchStates search_read_bwd(const Pattern &read,
                             const Pattern &kmer,
                             const KmerIndex &kmer_index,
                             const PRG_Info &prg_info);

SearchStates search_base_bwd(const Base &pattern_char,
                             const SearchStates &search_states,
                             const PRG_Info &prg_info);

AlleleId get_allele_id(const SA_Index &allele_marker_sa_index,
                       const PRG_Info &prg_info);

SA_Interval get_allele_marker_sa_interval(const Marker &site_marker_char,
                                          const PRG_Info &prg_info);

SearchStates process_markers_search_states(const SearchStates &search_states,
                                           const PRG_Info &prg_info);

SearchStates process_markers_search_state(const SearchState &search_state,
                                          const PRG_Info &prg_info);

#endif //GRAMTOOLS_SEARCH_HPP
