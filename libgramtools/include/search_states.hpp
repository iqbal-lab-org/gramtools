#include "prg.hpp"


#ifndef GRAMTOOLS_SEARCH_STATES_HPP
#define GRAMTOOLS_SEARCH_STATES_HPP

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

using SearchStates = std::list<SearchState>;

#endif //GRAMTOOLS_SEARCH_STATES_HPP
