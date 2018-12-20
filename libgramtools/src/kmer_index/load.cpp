#include <algorithm>
#include <thread>
#include <unordered_map>

#include "kmer_index/kmers.hpp"
#include "kmer_index/build.hpp"
#include "kmer_index/load.hpp"


using namespace gram;


Pattern gram::deserialize_next_kmer(const uint64_t &kmer_start_index,
                                    const sdsl::int_vector<3> &all_kmers,
                                    const uint32_t &kmers_size) {
    // TODO: implement as an iterator
    assert(kmer_start_index <= all_kmers.size() - kmers_size);
    Pattern kmer;
    kmer.reserve(kmers_size);
    for (uint64_t i = kmer_start_index; i < kmer_start_index + kmers_size; ++i)
        kmer.emplace_back(all_kmers[i]);
    return kmer;
}


IndexedKmerStats gram::deserialize_next_stats(const uint64_t &stats_index,
                                              const sdsl::int_vector<> &kmers_stats) {
    // TODO: implement as an iterator
    assert(stats_index < kmers_stats.size());
    IndexedKmerStats stats = {};

    stats.count_search_states = kmers_stats[stats_index];
    if (stats.count_search_states == 0)
        return stats;

    for (uint64_t i = stats_index + 1; i <= stats_index + stats.count_search_states; ++i)
        stats.path_lengths.push_back(kmers_stats[i]);
    return stats;
}


void pad_search_states(SearchStates &search_states,
                       const IndexedKmerStats &stats) {
    auto expected_count = stats.count_search_states - search_states.size();
    for (uint64_t i = search_states.size(); i < expected_count; ++i)
        search_states.emplace_back(SearchState{});
}


void handle_sa_interval(SearchStates &search_states,
                        uint64_t &sa_interval_index,
                        const sdsl::int_vector<> &sa_intervals) {
    for (auto &search_state: search_states) {
        search_state.sa_interval.first = sa_intervals[sa_interval_index];
        search_state.sa_interval.second = sa_intervals[sa_interval_index + 1];
        sa_interval_index += 2;
    }
}


void gram::parse_sa_intervals(KmerIndex &kmer_index,
                              const sdsl::int_vector<3> &all_kmers,
                              const sdsl::int_vector<> &kmers_stats,
                              const Parameters &parameters) {
    uint64_t sa_interval_index = 0;
    sdsl::int_vector<> sa_intervals;
    load_from_file(sa_intervals, parameters.sa_intervals_fpath);

    uint64_t stats_index = 0;
    uint64_t kmer_start_index = 0;

    while (kmer_start_index <= all_kmers.size() - parameters.kmers_size) {
        auto kmer = deserialize_next_kmer(kmer_start_index, all_kmers, parameters.kmers_size);
        kmer_start_index += parameters.kmers_size;

        auto stats = deserialize_next_stats(stats_index, kmers_stats);
        stats_index += stats.count_search_states + 1;

        auto &search_states = kmer_index[kmer];
        pad_search_states(search_states, stats);

        handle_sa_interval(search_states,
                           sa_interval_index,
                           sa_intervals);
    }
}


void handle_path_element(SearchStates &search_states,
                         uint64_t &paths_index,
                         const sdsl::int_vector<> &paths,
                         const IndexedKmerStats &stats) {
    uint64_t i = 0;
    for (auto &search_state: search_states) {
        const auto &path_length = stats.path_lengths[i++];

        for (uint64_t j = 0; j < path_length; ++j) {
            Marker marker = paths[paths_index];
            AlleleId allele_id = paths[paths_index + 1];
            paths_index += 2;

            VariantLocus site = {marker, allele_id};
            search_state.variant_site_path.emplace_back(site);
        }
    }
}


void gram::parse_paths(KmerIndex &kmer_index,
                       const sdsl::int_vector<3> &all_kmers,
                       const sdsl::int_vector<> &kmers_stats,
                       const Parameters &parameters) {
    uint64_t paths_index = 0;
    sdsl::int_vector<> paths;
    load_from_file(paths, parameters.paths_fpath);

    uint64_t stats_index = 0;
    uint64_t kmer_start_index = 0;

    while (kmer_start_index <= all_kmers.size() - parameters.kmers_size) {
        auto kmer = deserialize_next_kmer(kmer_start_index,
                                          all_kmers,
                                          parameters.kmers_size);
        kmer_start_index += parameters.kmers_size;

        auto stats = deserialize_next_stats(stats_index, kmers_stats);
        stats_index += stats.count_search_states + 1;

        auto &search_states = kmer_index[kmer];
        pad_search_states(search_states, stats);

        handle_path_element(search_states,
                            paths_index,
                            paths,
                            stats);
    }
}


KmerIndex gram::kmer_index::load(const Parameters &parameters) {
    KmerIndex kmer_index;

    sdsl::int_vector<3> all_kmers;
    load_from_file(all_kmers, parameters.kmers_fpath);

    sdsl::int_vector<> kmers_stats;
    load_from_file(kmers_stats, parameters.kmers_stats_fpath);

    parse_sa_intervals(kmer_index, all_kmers, kmers_stats, parameters);
    parse_paths(kmer_index, all_kmers, kmers_stats, parameters);
    return kmer_index;
}
