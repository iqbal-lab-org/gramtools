#include <sdsl/suffix_arrays.hpp>
#include "search.hpp"


SearchStates search_read_backwards(const Pattern &read,
                                   const Pattern &kmer,
                                   const KmerIndex &kmer_index,
                                   const PRG_Info &prg_info) {
    bool kmer_in_index = kmer_index.find(kmer) != kmer_index.end();
    if (not kmer_in_index)
        return SearchStates {};

    auto kmer_index_search_states = kmer_index.at(kmer);
    if (kmer_index_search_states.empty())
        return kmer_index_search_states;

    auto read_begin = read.rbegin();
    std::advance(read_begin, kmer.size());

    SearchStates new_search_states = kmer_index_search_states;

    for (auto it = read_begin; it != read.rend(); ++it) {
        const Base &pattern_char = *it;
        if (USE_SKIP_OPTIMIZATION) {
            set_state_skip_marker(new_search_states,
                                  prg_info);
        }
        new_search_states = process_read_char_search_states(pattern_char,
                                                            new_search_states,
                                                            prg_info);
        auto read_not_mapped = new_search_states.empty();
        if (read_not_mapped)
            break;
    }

    return new_search_states;
}


SearchState search_fm_index_base_backwards(const Base &pattern_char,
                                           const uint64_t char_first_sa_index,
                                           const SearchState &search_state,
                                           const PRG_Info &prg_info) {
    auto next_sa_interval = base_next_sa_interval(pattern_char,
                                                  char_first_sa_index,
                                                  search_state.sa_interval,
                                                  prg_info);
    auto valid_sa_interval = next_sa_interval.first - 1 != next_sa_interval.second;
    if (not valid_sa_interval) {
        SearchState new_search_state;
        new_search_state.invalid = true;
        return new_search_state;
    }

    auto new_search_state = search_state;
    new_search_state.sa_interval.first = next_sa_interval.first;
    new_search_state.sa_interval.second = next_sa_interval.second;
    return new_search_state;
}


SearchStates process_read_char_search_states(const Base &pattern_char,
                                             const SearchStates &old_search_states,
                                             const PRG_Info &prg_info) {
    auto post_markers_search_states = process_markers_search_states(old_search_states,
                                                                    prg_info);
    auto new_search_states = search_base_backwards(pattern_char,
                                                   post_markers_search_states,
                                                   prg_info);
    return new_search_states;
}


SearchState search_skipping_marker(const SearchState &search_state,
                                   const Base &pattern_char,
                                   const PRG_Info &prg_info) {
    SearchState new_search_state = search_state;

    auto previously_at_prg_edge = new_search_state.current_text_index == 0;
    if (previously_at_prg_edge) {
        new_search_state.invalid = true;
        return new_search_state;
    }

    new_search_state.current_text_index--;
    new_search_state.distance_to_next_marker--;

    auto searching_past_prg_start = new_search_state.current_text_index < 0;
    if (searching_past_prg_start) {
        new_search_state.invalid = true;
        return new_search_state;
    }

    auto search_is_valid = prg_info.encoded_prg[new_search_state.current_text_index]
                           == pattern_char;
    if (not search_is_valid) {
        new_search_state.invalid = true;
        return new_search_state;
    }

    auto stop_skipping_markers = new_search_state.distance_to_next_marker == 1;
    if (stop_skipping_markers) {
        new_search_state.skipping_to_marker = false;
        auto sa_index = prg_info.fm_index.isa[new_search_state.current_text_index];
        new_search_state.sa_interval.first = sa_index;
        new_search_state.sa_interval.second = sa_index;
        return new_search_state;
    }

    return new_search_state;
}


void set_state_skip_marker(SearchStates &search_states,
                           const PRG_Info &prg_info) {
    for (auto &search_state: search_states) {
        if (search_state.skipping_to_marker)
            continue;

        auto &sa_interval = search_state.sa_interval;
        auto state_unique_prg_match = sa_interval.first == sa_interval.second;
        if (not state_unique_prg_match)
            continue;

        auto &sa_index = sa_interval.first;
        auto state_text_index = prg_info.fm_index[sa_index];
        auto num_markers_before = prg_info.prg_markers_rank(state_text_index);

        if (num_markers_before == 0) {
            search_state.skipping_to_marker = true;
            search_state.distance_to_next_marker = state_text_index + 1;
            search_state.current_text_index = state_text_index;
            continue;
        }

        auto next_marker_index = prg_info.prg_markers_select(num_markers_before);
        auto distance_to_next_marker = state_text_index - next_marker_index;
        if (distance_to_next_marker > 1) {
            search_state.skipping_to_marker = true;
            search_state.distance_to_next_marker = distance_to_next_marker;
            search_state.current_text_index = state_text_index;
        }
    }
}


SA_Interval base_next_sa_interval(const Marker next_char,
                                  const SA_Index &next_char_first_sa_index,
                                  const SA_Interval &current_sa_interval,
                                  const PRG_Info &prg_info) {
    const auto &current_sa_start = current_sa_interval.first;
    const auto &current_sa_end = current_sa_interval.second;

    SA_Index sa_start_offset;
    if (current_sa_start <= 0)
        sa_start_offset = 0;
    else {
        if (next_char > 4)
            sa_start_offset = prg_info.fm_index.bwt.rank(current_sa_start, next_char);
        else {
            sa_start_offset = dna_bwt_rank(current_sa_start,
                                           next_char,
                                           prg_info);
        }
    }

    SA_Index sa_end_offset;
    if (next_char > 4)
        sa_end_offset = prg_info.fm_index.bwt.rank(current_sa_end + 1, next_char);
    else {
        sa_end_offset = dna_bwt_rank(current_sa_end + 1,
                                     next_char,
                                     prg_info);
    }

    auto new_start = next_char_first_sa_index + sa_start_offset;
    auto new_end = next_char_first_sa_index + sa_end_offset - 1;
    return SA_Interval {new_start, new_end};
}


void process_search_state_path_cache(SearchState &search_state) {
    if (not search_state.cache_populated)
        return;

    const auto &last_variant_site = search_state.variant_site_path.front();
    const auto cached_variant_site_already_recorded = search_state.cached_variant_site
                                                      == last_variant_site;
    if (not cached_variant_site_already_recorded)
        search_state.variant_site_path.push_front(search_state.cached_variant_site);

    search_state.cached_variant_site.first = 0;
    search_state.cached_variant_site.second = 0;
    search_state.cache_populated = false;
}


SearchStates search_base_backwards(const Base &pattern_char,
                                   const SearchStates &search_states,
                                   const PRG_Info &prg_info) {
    auto char_alphabet_rank = prg_info.fm_index.char2comp[pattern_char];
    auto char_first_sa_index = prg_info.fm_index.C[char_alphabet_rank];

    SearchStates new_search_states = {};

    for (const auto &search_state: search_states) {
        SearchState new_search_state;

        if (USE_SKIP_OPTIMIZATION and search_state.skipping_to_marker) {
            new_search_state = search_skipping_marker(search_state,
                                                      pattern_char,
                                                      prg_info);
        } else {
            new_search_state = search_fm_index_base_backwards(pattern_char,
                                                              char_first_sa_index,
                                                              search_state,
                                                              prg_info);
        }

        if (new_search_state.invalid)
            continue;
        process_search_state_path_cache(new_search_state);
        new_search_states.emplace_back(new_search_state);
    }

    return new_search_states;
}


SearchStates process_markers_search_states(const SearchStates &old_search_states,
                                           const PRG_Info &prg_info) {
    SearchStates new_search_states = old_search_states;
    SearchStates all_markers_new_search_states;
    for (const auto &search_state: old_search_states) {
        if (USE_SKIP_OPTIMIZATION and search_state.skipping_to_marker)
            continue;
        auto markers_search_states = process_markers_search_state(search_state, prg_info);
        all_markers_new_search_states.splice(all_markers_new_search_states.end(),
                                             markers_search_states);
    }
    new_search_states.splice(new_search_states.end(),
                             all_markers_new_search_states);
    return new_search_states;
}


struct SiteBoundaryMarkerInfo {
    bool is_start_boundary = false;
    SA_Interval sa_interval;
    Marker marker_char;
};


SiteBoundaryMarkerInfo site_boundary_marker_info(const Marker &marker_char,
                                                 const SA_Index &sa_right_of_marker,
                                                 const PRG_Info &prg_info) {
    auto alphabet_rank = prg_info.fm_index.char2comp[marker_char];
    auto first_sa_index = prg_info.fm_index.C[alphabet_rank];

    uint64_t marker_sa_index_offset;
    if (sa_right_of_marker <= 0)
        marker_sa_index_offset = 0;
    else
        marker_sa_index_offset = prg_info.fm_index.bwt.rank(sa_right_of_marker,
                                                            marker_char);
    auto marker_sa_index = first_sa_index + marker_sa_index_offset;
    const auto marker_text_idx = prg_info.fm_index[marker_sa_index];

    uint64_t other_marker_text_idx;
    if (marker_sa_index == first_sa_index)
        other_marker_text_idx = prg_info.fm_index[first_sa_index + 1];
    else
        other_marker_text_idx = prg_info.fm_index[first_sa_index];

    const bool marker_is_boundary_start = marker_text_idx <= other_marker_text_idx;
    return SiteBoundaryMarkerInfo {
            marker_is_boundary_start,
            SA_Interval {marker_sa_index, marker_sa_index},
            marker_char
    };
}


SA_Interval get_allele_marker_sa_interval(const Marker &site_marker_char,
                                          const PRG_Info &prg_info) {
    const auto allele_marker_char = site_marker_char + 1;
    const auto alphabet_rank = prg_info.fm_index.char2comp[allele_marker_char];
    const auto start_sa_index = prg_info.fm_index.C[alphabet_rank];

    const auto next_boundary_marker = allele_marker_char + 1;
    const auto max_char_in_alphabet = prg_info.fm_index.sigma - 1;
    const bool next_boundary_marker_valid = next_boundary_marker <= max_char_in_alphabet;

    SA_Index end_sa_index;
    if (next_boundary_marker_valid) {
        const auto next_boundary_marker_rank =
                prg_info.fm_index.char2comp[next_boundary_marker];
        const auto next_boundary_marker_start_sa_index =
                prg_info.fm_index.C[next_boundary_marker_rank];
        end_sa_index = next_boundary_marker_start_sa_index - 1;
    } else {
        end_sa_index = prg_info.fm_index.size() - 1;
    }
    return SA_Interval {start_sa_index, end_sa_index};
}


AlleleId get_allele_id(const SA_Index &allele_marker_sa_index,
                       const PRG_Info &prg_info) {
    auto internal_allele_text_index = prg_info.fm_index[allele_marker_sa_index] - 1;
    auto allele_id = (AlleleId) prg_info.allele_mask[internal_allele_text_index];
    return allele_id;
}


SearchStates get_allele_search_states(const Marker &site_boundary_marker,
                                      const SA_Interval &allele_marker_sa_interval,
                                      const SearchState &current_search_state,
                                      const PRG_Info &prg_info) {
    SearchStates search_states = {};

    const auto first_sa_interval_index = allele_marker_sa_interval.first;
    const auto last_sa_interval_index = allele_marker_sa_interval.second;

    for (auto allele_marker_sa_index = first_sa_interval_index;
         allele_marker_sa_index <= last_sa_interval_index;
         ++allele_marker_sa_index) {

        SearchState search_state = current_search_state;
        search_state.sa_interval.first = allele_marker_sa_index;
        search_state.sa_interval.second = allele_marker_sa_index;

        search_state.variant_site_state
                = SearchVariantSiteState::within_variant_site;
        search_state.cache_populated = true;

        auto allele_number = get_allele_id(allele_marker_sa_index, prg_info);
        search_state.cached_variant_site.first = site_boundary_marker;
        search_state.cached_variant_site.second = allele_number;

        search_states.emplace_back(search_state);
    }
    return search_states;
}


SearchState get_site_search_state(const AlleleId &final_allele_id,
                                  const SiteBoundaryMarkerInfo &boundary_marker_info,
                                  const SearchState &current_search_state,
                                  const PRG_Info &prg_info) {
    SearchState search_state = current_search_state;
    search_state.sa_interval.first = boundary_marker_info.sa_interval.first;
    search_state.sa_interval.second = boundary_marker_info.sa_interval.second;

    search_state.variant_site_state
            = SearchVariantSiteState::within_variant_site;

    search_state.cached_variant_site.first = boundary_marker_info.marker_char;
    search_state.cached_variant_site.second = final_allele_id;
    search_state.cache_populated = true;

    return search_state;
}


uint64_t get_number_of_alleles(const SA_Interval &allele_marker_sa_interval) {
    auto num_allele_markers = allele_marker_sa_interval.second
                              - allele_marker_sa_interval.first
                              + 1;
    auto num_alleles = num_allele_markers + 1;
    return num_alleles;
}


SearchStates entering_site_search_states(const SiteBoundaryMarkerInfo &boundary_marker_info,
                                         const SearchState &current_search_state,
                                         const PRG_Info &prg_info) {
    auto allele_marker_sa_interval =
            get_allele_marker_sa_interval(boundary_marker_info.marker_char,
                                          prg_info);
    auto new_search_states = get_allele_search_states(boundary_marker_info.marker_char,
                                                      allele_marker_sa_interval,
                                                      current_search_state,
                                                      prg_info);

    auto final_allele_id = get_number_of_alleles(allele_marker_sa_interval);
    auto site_search_state = get_site_search_state(final_allele_id,
                                                   boundary_marker_info,
                                                   current_search_state,
                                                   prg_info);
    new_search_states.emplace_back(site_search_state);
    return new_search_states;
}


SearchState exiting_site_search_state(const SiteBoundaryMarkerInfo &boundary_marker_info,
                                      const SearchState &current_search_state,
                                      const PRG_Info &prg_info) {
    SearchState new_search_state = current_search_state;
    bool read_started_in_allele = current_search_state.variant_site_state
                                  == SearchVariantSiteState::unknown;
    if (read_started_in_allele) {
        new_search_state.cache_populated = true;
        // allele ID is 1 because site boundary marker found
        // and exiting variant site with SearchVariantSiteState::unknown
        new_search_state.cached_variant_site.first = boundary_marker_info.marker_char;
        new_search_state.cached_variant_site.second = 1;
    }

    new_search_state.variant_site_state
            = SearchVariantSiteState::outside_variant_site;
    new_search_state.sa_interval.first = boundary_marker_info.sa_interval.first;
    new_search_state.sa_interval.second = boundary_marker_info.sa_interval.second;
    return new_search_state;
}


MarkersSearchResults left_markers_search(const SearchState &search_state,
                                         const PRG_Info &prg_info) {
    MarkersSearchResults markers_search_results;

    const auto &sa_interval = search_state.sa_interval;
    auto max_sa_index = sa_interval.second;
    auto sa_index = sa_interval.first;

    auto num_markers_before = prg_info.bwt_markers_rank(sa_index);
    uint64_t marker_count_offset = num_markers_before + 1;
    if (marker_count_offset > prg_info.markers_mask_count_set_bits)
        return markers_search_results;

    auto bwt_marker_index = prg_info.bwt_markers_select(marker_count_offset);

    while (bwt_marker_index <= max_sa_index) {
        auto marker = prg_info.fm_index.bwt[bwt_marker_index];
        auto search_result = std::make_pair(bwt_marker_index, marker);
        markers_search_results.emplace_back(search_result);

        ++marker_count_offset;
        if (marker_count_offset > prg_info.markers_mask_count_set_bits)
            return markers_search_results;
        bwt_marker_index = prg_info.bwt_markers_select(marker_count_offset);
    }

    return markers_search_results;
}


SearchStates process_boundary_marker(const Marker &marker_char,
                                     const SA_Index &sa_right_of_marker,
                                     const SearchState &current_search_state,
                                     const PRG_Info &prg_info) {
    auto boundary_marker_info = site_boundary_marker_info(marker_char,
                                                          sa_right_of_marker,
                                                          prg_info);

    bool entering_variant_site = not boundary_marker_info.is_start_boundary;
    if (entering_variant_site) {
        auto new_search_states = entering_site_search_states(boundary_marker_info,
                                                             current_search_state,
                                                             prg_info);
        return new_search_states;
    }

    bool exiting_variant_site = boundary_marker_info.is_start_boundary;
    if (exiting_variant_site) {
        auto new_search_state = exiting_site_search_state(boundary_marker_info,
                                                          current_search_state,
                                                          prg_info);
        return SearchStates {new_search_state};
    }
}


SearchState process_allele_marker(const Marker &allele_marker_char,
                                  const SA_Index &sa_right_of_marker,
                                  const SearchState &current_search_state,
                                  const PRG_Info &prg_info) {

    // end of allele found, skipping to variant site start boundary marker
    const Marker &boundary_marker_char = allele_marker_char - 1;

    auto alphabet_rank = prg_info.fm_index.char2comp[boundary_marker_char];
    auto first_sa_index = prg_info.fm_index.C[alphabet_rank];
    auto second_sa_index = first_sa_index + 1;

    SA_Index boundary_start_sa_index;
    bool boundary_start_is_first_sa = prg_info.fm_index[first_sa_index]
                                      < prg_info.fm_index[second_sa_index];
    if (boundary_start_is_first_sa)
        boundary_start_sa_index = first_sa_index;
    else
        boundary_start_sa_index = second_sa_index;

    auto new_search_state = current_search_state;
    new_search_state.sa_interval.first = boundary_start_sa_index;
    new_search_state.sa_interval.second = boundary_start_sa_index;
    new_search_state.variant_site_state = SearchVariantSiteState::outside_variant_site;

    auto internal_allele_text_index = prg_info.fm_index[sa_right_of_marker];
    auto allele_id = (AlleleId) prg_info.allele_mask[internal_allele_text_index];

    const auto &cached_variant_site = new_search_state.variant_site_path.front();
    bool read_started_within_allele = cached_variant_site.first != boundary_marker_char
                                      or cached_variant_site.second != allele_id;
    if (read_started_within_allele) {
        new_search_state.cached_variant_site.first = boundary_marker_char;
        new_search_state.cached_variant_site.second = allele_id;
        new_search_state.cache_populated = true;
    }
    return new_search_state;
}


SearchStates process_markers_search_state(const SearchState &current_search_state,
                                          const PRG_Info &prg_info) {
    const auto markers = left_markers_search(current_search_state,
                                             prg_info);
    if (markers.empty())
        return SearchStates {};

    SearchStates markers_search_states = {};

    for (const auto &marker: markers) {
        const auto &sa_right_of_marker = marker.first;
        const auto &marker_char = marker.second;

        const bool marker_is_site_boundary = marker_char % 2 == 1;
        if (marker_is_site_boundary) {
            auto new_search_states = process_boundary_marker(marker_char,
                                                             sa_right_of_marker,
                                                             current_search_state,
                                                             prg_info);
            markers_search_states.splice(markers_search_states.end(), new_search_states);
        } else {
            auto new_search_state = process_allele_marker(marker_char,
                                                          sa_right_of_marker,
                                                          current_search_state,
                                                          prg_info);
            markers_search_states.emplace_back(new_search_state);
        }
    }

    return markers_search_states;
}


std::string serialize_search_state(const SearchState &search_state) {
    std::stringstream ss;
    ss << "****** Search State ******" << std::endl;

    ss << "SA interval: ["
       << search_state.sa_interval.first
       << ", "
       << search_state.sa_interval.second
       << "]";
    ss << std::endl;

    if (not search_state.variant_site_path.empty()) {
        ss << "Variant site path [marker, allele id]: " << std::endl;
        for (const auto &variant_site: search_state.variant_site_path) {
            auto marker = variant_site.first;

            if (variant_site.second != 0) {
                const auto &allele_id = variant_site.second;
                ss << "[" << marker << ", " << allele_id << "]" << std::endl;
            }
        }
    }
    /*
    ss << "Variant site state: ";
    switch (search_state.variant_site_state) {
        case SearchVariantSiteState::within_variant_site:
            ss << "within_variant_site";
        case SearchVariantSiteState::outside_variant_site:
            ss << "outside_variant_site";
        case SearchVariantSiteState::unknown:
            ss << "unknown";
    }
    ss << std::endl;
     */

    switch (search_state.cache_populated) {
        case true:
            ss << "Cache populated: true" << std::endl;
            break;
        default:
            ss << "Cache populated: false" << std::endl;
            break;
    }
    if (search_state.cache_populated
        and search_state.cached_variant_site.first != 0) {
        const auto &marker = search_state.cached_variant_site.first;
        const auto &allele_id = search_state.cached_variant_site.second;
        ss << "Cached variant site: "
           << "[marker=" << marker
           << ", allele_id=" << allele_id
           << "]"
           << std::endl;
    }
    ss << "****** END Search State ******" << std::endl;
    return ss.str();
}


std::ostream &operator<<(std::ostream &os, const SearchState &search_state) {
    os << serialize_search_state(search_state);
    return os;
}
