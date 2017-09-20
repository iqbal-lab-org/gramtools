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

        bool process_markers = it != read_begin;
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


SearchStates search_char_bwd(const uint8_t pattern_char,
                             const SearchStates &search_states,
                             const PRG_Info &prg_info) {

    return SearchStates();


    /*
    const uint64_t num_begin =
        prg_info.fm_index.C[prg_info.fm_index.char2comp[marker_value - 1]];
    const auto site_boundary_lower = prg_info.fm_index[num_begin];
    const auto site_boundary_upper = prg_info.fm_index[num_begin + 1];

    if (site_boundary_lower < site_boundary_upper) {
        left = num_begin;
        right = num_begin + 1;
    } else {
        left = num_begin + 1;
        right = num_begin + 2;
    }
     */


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


SearchStates process_markers_search_state(const SearchState &search_state,
                                          const PRG_Info &prg_info) {

    const auto &sa_interval = search_state.sa_interval;
    const auto markers =
            prg_info.fm_index.wavelet_tree.range_search_2d(sa_interval.first,
                                                           sa_interval.second - 1,
                                                           5,
                                                           prg_info.max_alphabet_num).second;

    if (markers.empty())
        return SearchStates {search_state};

    SearchStates markers_search_states = {};
    uint64_t marker_boundary_char = 0;

    for (const auto &marker: markers) {
        const auto marker_index = marker.first;
        const auto marker_char = marker.second;

        const bool marker_is_boundary_char = marker_char % 2 == 1;
        if (marker_is_boundary_char)
            marker_boundary_char = marker_char;
        else
            marker_boundary_char = marker_char - 1;
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
