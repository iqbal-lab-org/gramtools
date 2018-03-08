#include <algorithm>
#include <thread>
#include <unordered_map>

#include "search/search.hpp"
#include "kmer_index/load.hpp"
#include "kmer_index/dump.hpp"


KmerIndexStats calculate_stats(const KmerIndex &kmer_index) {
    KmerIndexStats stats = {};
    stats.count_kmers = kmer_index.size();

    for (const auto &entry: kmer_index) {
        // memory elements for recording path lengths of each search state
        const auto &search_states = entry.second;
        stats.count_search_states += search_states.size();

        for (const auto &search_state: search_states)
            stats.count_total_path_elements += search_state.variant_site_path.size() * 2;
    }
    return stats;
}

sdsl::int_vector<3> dump_kmers(const KmerIndex &kmer_index,
                               const Parameters &parameters) {
    sdsl::int_vector<3> all_kmers(kmer_index.size() * parameters.kmers_size);
    uint64_t i = 0;

    for (const auto &entry: kmer_index) {
        const auto &kmer = entry.first;
        for (const auto &base: kmer) {
            assert(base >= 1 and base <= 4);
            all_kmers[i++] = base;
        }
    }

    store_to_file(all_kmers, parameters.kmers_fpath);
    return all_kmers;
}

void dump_kmers_stats(const KmerIndexStats &stats,
                      const sdsl::int_vector<3> &all_kmers,
                      const KmerIndex &kmer_index,
                      const Parameters &parameters) {
    // each kmer: number of search states, path length, path length...
    const auto &count_distinct_paths = stats.count_search_states;
    uint64_t count_memory_elements = stats.count_kmers + count_distinct_paths;
    sdsl::int_vector<> kmers_stats(count_memory_elements, 1, 16);
    uint64_t i = 0;

    uint64_t kmer_start_index = 0;
    while (kmer_start_index <= all_kmers.size() - parameters.kmers_size) {
        auto kmer = deserialize_next_kmer(kmer_start_index,
                                          all_kmers,
                                          parameters.kmers_size);
        kmer_start_index += kmer.size();

        const auto &search_states = kmer_index.at(kmer);
        kmers_stats[i++] = search_states.size();
        for (const auto &search_state: search_states)
            kmers_stats[i++] = search_state.variant_site_path.size();
    }

    sdsl::util::bit_compress(kmers_stats);
    store_to_file(kmers_stats, parameters.kmers_stats_fpath);
}

void dump_sa_intervals(const KmerIndexStats &stats,
                       const sdsl::int_vector<3> &all_kmers,
                       const KmerIndex &kmer_index,
                       const Parameters &parameters) {
    sdsl::int_vector<> sa_intervals(stats.count_search_states * 2, 0, 32);
    uint64_t i = 0;

    uint64_t kmer_start_index = 0;
    while (kmer_start_index <= all_kmers.size() - parameters.kmers_size) {
        auto kmer = deserialize_next_kmer(kmer_start_index,
                                          all_kmers,
                                          parameters.kmers_size);
        kmer_start_index += kmer.size();

        const auto &search_states = kmer_index.at(kmer);
        for (const auto &search_state: search_states) {
            sa_intervals[i++] = search_state.sa_interval.first;
            sa_intervals[i++] = search_state.sa_interval.second;
        }
    }

    sdsl::util::bit_compress(sa_intervals);
    store_to_file(sa_intervals, parameters.sa_intervals_fpath);
}

void dump_paths(const KmerIndexStats &stats,
                const sdsl::int_vector<3> &all_kmers,
                const KmerIndex &kmer_index,
                const Parameters &parameters) {
    sdsl::int_vector<> paths(stats.count_total_path_elements, 0, 32);
    uint64_t i = 0;

    uint64_t kmer_start_index = 0;
    while (kmer_start_index <= all_kmers.size() - parameters.kmers_size) {
        auto kmer = deserialize_next_kmer(kmer_start_index,
                                          all_kmers,
                                          parameters.kmers_size);
        kmer_start_index += kmer.size();

        const auto &search_states = kmer_index.at(kmer);
        for (const auto &search_state: search_states) {
            for (const auto &path_element: search_state.variant_site_path) {
                paths[i++] = path_element.first;
                paths[i++] = path_element.second;
            }
        }
    }

    sdsl::util::bit_compress(paths);
    store_to_file(paths, parameters.paths_fpath);
}

void kmer_index::dump(const KmerIndex &kmer_index,
                      const Parameters &parameters) {
    sdsl::int_vector<3> all_kmers = dump_kmers(kmer_index, parameters);
    auto stats = calculate_stats(kmer_index);
    dump_kmers_stats(stats, all_kmers, kmer_index, parameters);
    dump_sa_intervals(stats, all_kmers, kmer_index, parameters);
    dump_paths(stats, all_kmers, kmer_index, parameters);
}