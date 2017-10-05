#include "kmers.hpp"


#ifndef GRAMTOOLS_BIDIR_SEARCH_BWD_V2_HPP
#define GRAMTOOLS_BIDIR_SEARCH_BWD_V2_HPP


enum class SearchVariantSiteState {
    within_variant_site,
    outside_variant_site,
    unknown
};

struct SearchState {
    SA_Interval sa_interval;
    VariantSitePath variant_site_path;

    SearchVariantSiteState variant_site_state = SearchVariantSiteState::unknown;
    bool variant_site_recorded = false;
    VariantSite last_variant_site;

    bool operator==(const SearchState &other) const {
        return this->sa_interval == other.sa_interval
               and this->variant_site_path == other.variant_site_path
               and this->variant_site_state == other.variant_site_state
               and this->variant_site_recorded == other.variant_site_recorded
               and this->last_variant_site == other.last_variant_site;
    };
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

uint64_t get_allele_id(const uint64_t allele_marker_sa_index,
                       const PRG_Info &prg_info);

SA_Interval get_allele_marker_sa_interval(const uint64_t site_marker_char,
                                          const PRG_Info &prg_info);

SearchStates process_markers_search_states(const SearchStates &search_states,
                                           const PRG_Info &prg_info);

SearchStates process_markers_search_state(const SearchState &search_state,
                                          const PRG_Info &prg_info);


#endif //GRAMTOOLS_BIDIR_SEARCH_BWD_V2_HPP
