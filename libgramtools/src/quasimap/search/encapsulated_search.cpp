#include "quasimap/search/encapsulated_search.hpp"

/**
 * A caching object used to temporarily store a single search state
 * @see handle_allele_encapsulated_state()
 */
class SearchStateCache {
public:
    SearchState search_state = {};
    bool empty = true;

    void set(const SearchState &search_state) {
        this->search_state = search_state;
        this->empty = false;
    }

    void flush(SearchStates &search_states) {
        if (this->empty)
            return;
        search_states.emplace_back(this->search_state);
        this->empty = true;
    }

    void update_sa_interval_max(const SA_Index &new_sa_interval_max) {
        assert(not this->empty);
        assert(this->search_state.sa_interval.second + 1 == new_sa_interval_max);
        this->search_state.sa_interval.second = new_sa_interval_max;
    }
};


SearchStates gram::handle_allele_encapsulated_state(const SearchState &search_state,
                                                    const PRG_Info &prg_info) {
    bool has_path = not search_state.traversed_path.empty();
    assert(not has_path);

    SearchStates new_search_states = {};
    SearchStateCache cache;

    for (uint64_t sa_index = search_state.sa_interval.first;
         sa_index <= search_state.sa_interval.second;
         ++sa_index) {

        auto prg_index = prg_info.fm_index[sa_index];
        auto site_marker = prg_info.sites_mask[prg_index];
        auto allele_id = prg_info.allele_mask[prg_index];

        bool within_site = site_marker != 0;
        if (not within_site) {
            cache.flush(new_search_states);
            cache.set(SearchState{
                    SA_Interval{sa_index, sa_index},
                    VariantSitePath{},
                    VariantSitePath{},
                    SearchVariantSiteState::outside_variant_site
            });
            cache.flush(new_search_states);
            continue;
        }

        //  else: read is completely encapsulated within allele

        if (cache.empty) {
            cache.set(SearchState{
                    SA_Interval{sa_index, sa_index},
                    VariantSitePath{
                            VariantLocus{site_marker, allele_id}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            });
            continue;
        }

        // cache is not empty: if the current SearchState goes through same site and allele,
        // lengthen the SA interval of cached SearchState
        VariantSitePath current_path = {VariantLocus{site_marker, allele_id}};
        bool cache_has_same_path = current_path == cache.search_state.traversed_path;
        if (cache_has_same_path) {
            cache.update_sa_interval_max(sa_index);
            continue;
        } else {
            cache.flush(new_search_states);
            cache.set(SearchState{
                    SA_Interval{sa_index, sa_index},
                    current_path,
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            });
        }
    }
    cache.flush(new_search_states);
    return new_search_states;
}

SearchStates gram::handle_allele_encapsulated_states(const SearchStates &search_states,
                                                     const PRG_Info &prg_info) {
    SearchStates new_search_states = {};

    for (const auto &search_state: search_states) {
        bool has_path = not search_state.traversed_path.empty();
        if (has_path) {
            new_search_states.emplace_back(search_state);
            continue;
        }

        SearchStates split_search_states = handle_allele_encapsulated_state(search_state,
                                                                            prg_info);
        for (const auto &split_search_state: split_search_states)
            new_search_states.emplace_back(split_search_state);
    }
    return new_search_states;
}
