#include "utils.hpp"
#include "prg.hpp"
#include "kmer_index.hpp"


#ifndef GRAMTOOLS_SEARCH_HPP
#define GRAMTOOLS_SEARCH_HPP

enum class SearchVariantSiteState {
    within_variant_site,
    outside_variant_site,
    unknown
};

struct SearchState {
    SA_Interval sa_interval;
    VariantSitePath variant_site_path;

    SearchVariantSiteState variant_site_state = SearchVariantSiteState::unknown;

    bool cache_populated = false;
    VariantSite cached_variant_site;

    bool operator== (const SearchState &other) const {
        return this->sa_interval == other.sa_interval
               and this->variant_site_path == other.variant_site_path
               and this->variant_site_state == other.variant_site_state
               and this->cache_populated == other.cache_populated
               and this->cached_variant_site == other.cached_variant_site;
    };
};

std::string serialize_search_state(const SearchState &search_state);

std::ostream &operator<< (std::ostream &os, const SearchState &search_state);

using SearchStates = std::list<SearchState>;

SearchStates get_initial_search_states(const Pattern &kmer,
                                       const KmerIndex &kmer_index);

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
