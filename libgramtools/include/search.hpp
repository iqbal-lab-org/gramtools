#include "kmers.hpp"


#ifndef GRAMTOOLS_BIDIR_SEARCH_BWD_V2_HPP
#define GRAMTOOLS_BIDIR_SEARCH_BWD_V2_HPP


struct SearchState {
    SA_Interval sa_interval;
    Site site;
};

using SearchStates = std::list<SearchState>;
using Pattern = std::vector<uint8_t>;

SearchStates initial_search_states(const Pattern &kmer,
                                   const KmerIndex &kmer_index);

SearchStates search_read_bwd(const Pattern &read,
                             const Pattern &kmer,
                             const KmerIndex &kmer_index,
                             const PRG_Info &prg_info);

SearchStates search_char_bwd(const uint8_t pattern_char,
                             const SearchStates &search_states,
                             const PRG_Info &prg_info);

SearchStates process_markers_search_states(const SearchStates &search_states,
                                           const PRG_Info &prg_info);

SearchStates process_markers_search_state(const SearchState &search_state,
                                          const PRG_Info &prg_info);


#endif //GRAMTOOLS_BIDIR_SEARCH_BWD_V2_HPP
