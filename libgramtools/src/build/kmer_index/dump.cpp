#include <algorithm>
#include <thread>
#include <unordered_map>

#include "build/kmer_index/load.hpp"
#include "build/kmer_index/dump.hpp"


using namespace gram;


KmerIndexStats gram::calculate_stats(const KmerIndex &kmer_index) {
    KmerIndexStats stats = {};
    stats.count_kmers = kmer_index.size();

    for (const auto &entry: kmer_index) {
        // memory elements for recording path lengths of each search state
        const auto &search_states = entry.second;
        stats.count_search_states += search_states.size();

        for (const auto &search_state: search_states)
            stats.count_total_path_elements +=
                    2 * search_state.traversing_path.size() +
                    2 * search_state.traversed_path.size();
    }
    return stats;
}


sdsl::int_vector<3> gram::dump_kmers(const KmerIndex &kmer_index,
                                     const BuildParams &parameters) {
    //Constructor parameter passed: total number of bases to store.
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


void gram::dump_kmers_stats(const KmerIndexStats &stats,
                            const sdsl::int_vector<3> &all_kmers,
                            const KmerIndex &kmer_index,
                            const BuildParams &parameters) {
    // Makes room for storing for each kmer: number of search states, path length, path length...
    const auto &count_distinct_paths = stats.count_search_states;
    // For each kmer, we will write the number of search states it has,
    // and for each of those, we will write the number of `VariantLocus` (loci)
    uint64_t count_memory_elements = stats.count_kmers + count_distinct_paths;
    sdsl::int_vector<> kmers_stats(count_memory_elements, 0, 32);
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
            kmers_stats[i++] =
                    search_state.traversing_path.size() +
                    search_state.traversed_path.size();
    }

    sdsl::util::bit_compress(kmers_stats);
    store_to_file(kmers_stats, parameters.kmers_stats_fpath);
}


void gram::dump_sa_intervals(const KmerIndexStats &stats,
                             const sdsl::int_vector<3> &all_kmers,
                             const KmerIndex &kmer_index,
                             const BuildParams &parameters) {
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


void gram::dump_paths(const KmerIndexStats &stats,
                      const sdsl::int_vector<3> &all_kmers,
                      const KmerIndex &kmer_index,
                      const BuildParams &parameters) {
    // Sdsl stores unsigned integer vectors, so scale Allele IDs up to >0 if needed
    AlleleId increment{0};
    if (ALLELE_UNKNOWN < 0) increment = std::abs(ALLELE_UNKNOWN);

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
            for (const auto &path_element: search_state.traversed_path) {
                paths[i++] = path_element.first;
                paths[i++] = path_element.second + increment;
            }
            for (const auto &path_element: search_state.traversing_path) {
                assert(path_element.second == ALLELE_UNKNOWN);
                paths[i++] = path_element.first;
                paths[i++] = path_element.second + increment;
            }
        }
    }

    sdsl::util::bit_compress(paths);
    store_to_file(paths, parameters.paths_fpath);
}


void gram::kmer_index::dump(const KmerIndex &kmer_index,
                            const BuildParams &parameters) {
    sdsl::int_vector<3> all_kmers = dump_kmers(kmer_index, parameters);
    auto stats = calculate_stats(kmer_index);
    dump_kmers_stats(stats, all_kmers, kmer_index, parameters);
    dump_sa_intervals(stats, all_kmers, kmer_index, parameters);
    dump_paths(stats, all_kmers, kmer_index, parameters);
}