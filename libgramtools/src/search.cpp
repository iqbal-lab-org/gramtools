#include <sdsl/suffix_arrays.hpp>

#include "kmers.hpp"
#include "map.hpp"
#include "search.hpp"


SearchStates search_read_bwd(const Pattern &read,
                             const Pattern &kmer,
                             const KmerIndex &kmer_index,
                             const PRG_Info &prg_info) {

    auto search_states = initial_search_states(kmer, kmer_index);
    if (search_states.empty())
        return search_states;

    auto read_begin = read.rbegin();
    std::advance(read_begin, kmer.size());

    for (auto it = read_begin; it != read.rend(); ++it) {
        SearchStates new_search_states;

        const bool process_markers = it != read_begin;
        if (process_markers)
            new_search_states = process_markers_search_states(search_states,
                                                              prg_info);
        else
            new_search_states = search_states;

        const uint8_t read_char = *it;
        search_states = search_char_bwd(read_char,
                                        new_search_states,
                                        prg_info);
        if (search_states.empty())
            break;
    }

    return search_states;
}


SA_Interval get_next_sa_interval(const uint64_t current_char,
                                 const uint64_t current_char_first_sa_index,
                                 const SA_Interval &current_sa_interval,
                                 const PRG_Info &prg_info) {
    const auto &current_sa_start = current_sa_interval.first;
    const auto &current_sa_end = current_sa_interval.second;
    assert(current_sa_end <= prg_info.fm_index.size() - 1);

    uint64_t sa_start_offset;
    if (current_sa_start <= 0)
        sa_start_offset = 0;
    else
        sa_start_offset = prg_info.fm_index.bwt.rank(current_sa_start, current_char);

    uint64_t sa_end_offset = prg_info.fm_index.bwt.rank(current_sa_end + 1, current_char);

    auto new_start = current_char_first_sa_index + sa_start_offset;
    auto new_end = current_char_first_sa_index + sa_end_offset - 1;
    return SA_Interval {new_start, new_end};
}


SearchStates search_char_bwd(const uint8_t pattern_char,
                             const SearchStates &search_states,
                             const PRG_Info &prg_info) {
    auto char_alphabet_rank = prg_info.fm_index.char2comp[pattern_char];
    auto char_first_sa_index = prg_info.fm_index.C[char_alphabet_rank];

    SearchStates new_search_states = {};

    for (const auto &search_state: search_states) {
        auto next_sa_interval = get_next_sa_interval(pattern_char,
                                                     char_first_sa_index,
                                                     search_state.sa_interval,
                                                     prg_info);
        auto valid_sa_interval = next_sa_interval.first - 1 != next_sa_interval.second;
        if (!valid_sa_interval)
            continue;

        // TODO: base next_search_state on search_state

        SearchState next_search_state = {
                next_sa_interval,
                search_state.variant_site_path
        };
        new_search_states.emplace_back(next_search_state);
    }

    return new_search_states;
}


SearchStates process_markers_search_states(const SearchStates &search_states,
                                           const PRG_Info &prg_info) {
    SearchStates new_search_states;
    for (const auto &search_state: search_states) {
        auto markers_search_states = process_markers_search_state(search_state, prg_info);
        new_search_states.splice(new_search_states.end(), markers_search_states);
    }
    return new_search_states;
}


struct SiteBoundaryMarkerInfo {
    bool is_start_boundary = false;
    SA_Interval sa_interval;
    uint64_t marker_char;
};


SiteBoundaryMarkerInfo site_boundary_marker_info(const uint64_t &marker_char,
                                                 const uint64_t &sa_preceding_marker_index,
                                                 const PRG_Info &prg_info) {
    auto alphabet_rank = prg_info.fm_index.char2comp[marker_char];
    auto first_sa_index = prg_info.fm_index.C[alphabet_rank];

    uint64_t marker_sa_index_offset;
    if (sa_preceding_marker_index <= 0)
        marker_sa_index_offset = 0;
    else
        marker_sa_index_offset = prg_info.fm_index.bwt.rank(sa_preceding_marker_index,
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


SA_Interval get_allele_marker_sa_interval(const uint64_t site_marker_char,
                                          const PRG_Info &prg_info) {
    const auto allele_marker_char = site_marker_char + 1;
    const auto alphabet_rank = prg_info.fm_index.char2comp[allele_marker_char];

    const auto start_sa_index = prg_info.fm_index.C[alphabet_rank];
    uint64_t end_sa_index = start_sa_index + 1;
    return SA_Interval {start_sa_index, end_sa_index};
}


uint64_t get_allele_id(const uint64_t allele_marker_sa_index,
                       const PRG_Info &prg_info) {
    auto internal_allele_text_index = prg_info.fm_index[allele_marker_sa_index] - 1;
    auto allele_id = (uint64_t) prg_info.allele_mask[internal_allele_text_index];
    return allele_id;
}


SearchStates get_allele_search_states(const uint64_t site_boundary_marker,
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
        search_state.cached_variant_site.second = Allele {allele_number};

        search_states.emplace_back(search_state);
    }
    return search_states;
}


SearchState get_site_search_state(const uint64_t final_allele_id,
                                  const SiteBoundaryMarkerInfo &boundary_marker_info,
                                  const SearchState &current_search_state,
                                  const PRG_Info &prg_info) {
    SearchState search_state = current_search_state;
    search_state.sa_interval.first = boundary_marker_info.sa_interval.first;
    search_state.sa_interval.second = boundary_marker_info.sa_interval.second;

    search_state.variant_site_state
            = SearchVariantSiteState::within_variant_site;

    search_state.cached_variant_site.first = boundary_marker_info.marker_char;
    search_state.cached_variant_site.second = Allele {final_allele_id};
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
        // allele 1 because site boundary marker found
        // and exiting variant site with SearchVariantSiteState::unknown
        new_search_state.cached_variant_site = {boundary_marker_info.marker_char, Allele {1}};
    }

    new_search_state.variant_site_state
            = SearchVariantSiteState::outside_variant_site;
    new_search_state.sa_interval.first = boundary_marker_info.sa_interval.first;
    new_search_state.sa_interval.second = boundary_marker_info.sa_interval.second;
    return new_search_state;
}


using MarkerIndexPrecedingSA = uint64_t;
using MarkerChar = uint64_t;
using MarkersSearchResults = std::vector<std::pair<MarkerIndexPrecedingSA, MarkerChar>>;


MarkersSearchResults left_markers_search(const SearchState &search_state,
                                         const PRG_Info &prg_info) {
    const auto &sa_interval = search_state.sa_interval;
    const auto &wt = prg_info.fm_index.wavelet_tree;
    const auto wt_search_result = wt.range_search_2d(sa_interval.first, sa_interval.second,
                                                     5, prg_info.max_alphabet_num);
    const auto &markers = wt_search_result.second;
    return markers;
}


SearchStates process_boundary_marker(const uint64_t marker_char,
                                     const uint64_t sa_preceding_marker_index,
                                     const SearchState &current_search_state,
                                     const PRG_Info &prg_info) {
    auto boundary_marker_info = site_boundary_marker_info(marker_char,
                                                          sa_preceding_marker_index,
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


SearchState process_allele_marker(const uint64_t allele_marker_char,
                                  const uint64_t sa_preceding_marker_index,
                                  const SearchState &current_search_state,
                                  const PRG_Info &prg_info) {

    // end of allele found, skipping to variant site start boundary marker
    const uint64_t boundary_marker_char = allele_marker_char - 1;

    auto alphabet_rank = prg_info.fm_index.char2comp[boundary_marker_char];
    auto first_sa_index = prg_info.fm_index.C[alphabet_rank];
    auto second_sa_index = first_sa_index + 1;

    uint64_t boundary_start_sa_index;
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

    auto internal_allele_text_index = prg_info.fm_index[sa_preceding_marker_index];
    auto allele_id = (uint64_t) prg_info.allele_mask[internal_allele_text_index];

    const auto &cached_variant_site = new_search_state.variant_site_path.front();
    bool read_started_within_allele = cached_variant_site.first != boundary_marker_char
                                      or cached_variant_site.second != Allele {allele_id};
    if (read_started_within_allele) {
        new_search_state.cached_variant_site = {boundary_marker_char, Allele {allele_id}};
        new_search_state.cache_populated = true;
    }
    return new_search_state;
}


SearchStates process_markers_search_state(const SearchState &current_search_state,
                                          const PRG_Info &prg_info) {
    const auto markers = left_markers_search(current_search_state,
                                             prg_info);

    if (markers.empty())
        return SearchStates {current_search_state};

    SearchStates markers_search_states = {};

    for (const auto &marker: markers) {
        const auto &sa_preceding_marker_index = marker.first;
        const auto &marker_char = marker.second;

        const bool marker_is_site_boundary = marker_char % 2 == 1;
        if (marker_is_site_boundary) {
            auto new_search_states = process_boundary_marker(marker_char,
                                                             sa_preceding_marker_index,
                                                             current_search_state,
                                                             prg_info);
            markers_search_states.splice(markers_search_states.end(), new_search_states);
        } else {
            auto new_search_state = process_allele_marker(marker_char,
                                                          sa_preceding_marker_index,
                                                          current_search_state,
                                                          prg_info);
            markers_search_states.emplace_back(new_search_state);
        }
    }

    return markers_search_states;
}


SearchStates initial_search_states(const Pattern &kmer,
                                   const KmerIndex &kmer_index) {

    SearchStates search_states = {};
    const bool kmer_not_in_index = kmer_index.sa_intervals_map.find(kmer)
                                   == kmer_index.sa_intervals_map.end();
    if (kmer_not_in_index)
        return search_states;

    const auto &sa_intervals = kmer_index.sa_intervals_map.at(kmer);
    const auto &sites = kmer_index.sites_map.at(kmer);

    auto sa_intervals_it = sa_intervals.begin();
    auto sites_it = sites.begin();
    while (sa_intervals_it != sa_intervals.end()) {
        auto sa_interval = *sa_intervals_it;
        auto site = *sites_it;
        SearchState search_state = {
                sa_interval,
                site
        };
        search_states.emplace_back(search_state);

        ++sa_intervals_it;
        ++sites_it;
    }

    return search_states;
}
