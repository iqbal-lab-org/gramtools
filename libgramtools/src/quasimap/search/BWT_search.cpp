#include <sdsl/suffix_arrays.hpp>
#include "quasimap/search/BWT_search.hpp"


using namespace gram;

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

void gram::set_allele_ids(SearchStates &search_states,
                            const PRG_Info &prg_info) {

    for (auto &search_state: search_states) {
        // If the variant site path is empty, we cannot have an unknown allele id.
        if (!search_state.traversing_path.empty()) {

            auto last = search_state.traversing_path.back();
            search_state.traversing_path.pop_back();
            assert(search_state.traversing_path.size() == 0);

            for (int SA_pos = search_state.sa_interval.first; SA_pos <= search_state.sa_interval.second; SA_pos++) {
                auto new_search_state = search_state;
                // Find out allele id
                auto text_pos = prg_info.fm_index[SA_pos];
                auto allele_id = prg_info.allele_mask[text_pos];
                last.second = allele_id;

                new_search_state.traversed_path.emplace_back(last);
                new_search_state.sa_interval = SA_Interval{SA_pos, SA_pos};

                if (SA_pos != search_state.sa_interval.second){ // Case: add to the set of search states
                    search_states.emplace_back(new_search_state);
                } 
                else search_state = new_search_state; // Case: end of the iteration; modify the search state in place.
            }
        }
    }

}


/**
 * Backward search followed by check whether the extended searched pattern maps somewhere in the prg.
 */
SearchState search_fm_index_base_backwards(const int_Base &pattern_char,
                                           const uint64_t char_first_sa_index,
                                           const SearchState &search_state,
                                           const PRG_Info &prg_info) {
    auto next_sa_interval = base_next_sa_interval(pattern_char,
                                                  char_first_sa_index,
                                                  search_state.sa_interval,
                                                  prg_info);
    //  An 'invalid' SA interval (i,j) is defined by i-1=j, which occurs when the read no longer maps anywhere in the prg.
    auto valid_sa_interval = next_sa_interval.first - 1 != next_sa_interval.second;
    if (not valid_sa_interval) { // Create an empty, invalid search state.
        SearchState new_search_state;
        new_search_state.invalid = true;
        return new_search_state;
    }

    auto new_search_state = search_state;
    new_search_state.sa_interval.first = next_sa_interval.first;
    new_search_state.sa_interval.second = next_sa_interval.second;
    return new_search_state;
}


SA_Interval gram::base_next_sa_interval(const Marker &next_char,
                                        const SA_Index &next_char_first_sa_index,
                                        const SA_Interval &current_sa_interval,
                                        const PRG_Info &prg_info) {
    const auto &current_sa_start = current_sa_interval.first;
    const auto &current_sa_end = current_sa_interval.second;

    SA_Index sa_start_offset;
    if (current_sa_start <= 0)
        sa_start_offset = 0;
    else {
        //  TODO: Consider deleting this if-clause, next_char should never be > 4, it probably never runs
        if (next_char > 4)
            sa_start_offset = prg_info.fm_index.bwt.rank(current_sa_start, next_char);
        else {
            sa_start_offset = dna_bwt_rank(current_sa_start,
                                           next_char,
                                           prg_info);
        }
    }

    SA_Index sa_end_offset;
    //  TODO: Consider deleting this if-clause, next_char should never be > 4, it probably never runs
    if (next_char > 4)
        sa_end_offset = prg_info.fm_index.bwt.rank(current_sa_end + 1, next_char);
    else {
        sa_end_offset = dna_bwt_rank(current_sa_end + 1,
                                     next_char,
                                     prg_info);
    }

    auto new_start = next_char_first_sa_index + sa_start_offset;
    auto new_end = next_char_first_sa_index + sa_end_offset - 1;
    return SA_Interval{new_start, new_end};
}

SearchStates gram::search_base_backwards(const int_Base &pattern_char,
                                         const SearchStates &search_states,
                                         const PRG_Info &prg_info) {
    // Compute the first occurrence of `pattern_char` in the suffix array. Necessary for backward search.
    auto char_alphabet_rank = prg_info.fm_index.char2comp[pattern_char];
    auto char_first_sa_index = prg_info.fm_index.C[char_alphabet_rank];

    SearchStates new_search_states = {};

    for (const auto &search_state: search_states) {
        SearchState new_search_state = search_fm_index_base_backwards(pattern_char,
                                                                      char_first_sa_index,
                                                                      search_state,
                                                                      prg_info);
        if (new_search_state.invalid)
            continue;
        new_search_states.emplace_back(new_search_state);
    }

    return new_search_states;
}




std::string gram::serialize_search_state(const SearchState &search_state) {
    std::stringstream ss;
    ss << "****** Search State ******" << std::endl;

    ss << "SA interval: ["
       << search_state.sa_interval.first
       << ", "
       << search_state.sa_interval.second
       << "]";
    ss << std::endl;

    if (not search_state.traversed_path.empty()) {
        ss << "Variant site path [marker, allele id]: " << std::endl;
        for (const auto &variant_site: search_state.traversed_path) {
            auto marker = variant_site.first;

            if (variant_site.second != 0) {
                const auto &allele_id = variant_site.second;
                ss << "[" << marker << ", " << allele_id << "]" << std::endl;
            }
        }
    }
    ss << "****** END Search State ******" << std::endl;
    return ss.str();
}


std::ostream &gram::operator<<(std::ostream &os, const SearchState &search_state) {
    os << serialize_search_state(search_state);
    return os;
}
