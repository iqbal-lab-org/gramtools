#include <algorithm>
#include <thread>
#include <unordered_map>
#include <kmer_index/kmers.hpp>

#include "kmer_index/kmer_index.hpp"
#include "search/search.hpp"


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
        for (const auto &base: kmer)
            all_kmers[i++] = base;
    }

    sdsl::store_to_file(all_kmers, parameters.kmers_fpath);
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
    sdsl::store_to_file(kmers_stats, parameters.kmers_stats_fpath);
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
    sdsl::store_to_file(sa_intervals, parameters.sa_intervals_fpath);
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
    sdsl::store_to_file(paths, parameters.paths_fpath);
}


void dump_kmer_index(const KmerIndex &kmer_index,
                     const Parameters &parameters) {
    sdsl::int_vector<3> all_kmers = dump_kmers(kmer_index, parameters);
    auto stats = calculate_stats(kmer_index);
    dump_kmers_stats(stats, all_kmers, kmer_index, parameters);
    dump_sa_intervals(stats, all_kmers, kmer_index, parameters);
    dump_paths(stats, all_kmers, kmer_index, parameters);
}


CacheElement get_next_cache_element(const Base &base,
                                    const bool kmer_base_is_last,
                                    const CacheElement &last_cache_element,
                                    const PRG_Info &prg_info) {
    const auto &old_search_states = last_cache_element.search_states;
    SearchStates new_search_states;
    if (not kmer_base_is_last) {
        new_search_states = process_markers_search_states(old_search_states,
                                                          prg_info);
    } else {
        new_search_states = old_search_states;
    }
    new_search_states = search_base_backwards(base,
                                              new_search_states,
                                              prg_info);
    return CacheElement {
            new_search_states,
            base
    };
}


CacheElement get_initial_cache_element(const Base &base,
                                       const PRG_Info &prg_info) {
    SearchState search_state = {
            SA_Interval {0, prg_info.fm_index.size() - 1}
    };
    SearchStates search_states = {search_state};
    CacheElement initial_cache_element = {search_states};

    bool kmer_base_is_last = true;
    const auto &cache_element = get_next_cache_element(base,
                                                       kmer_base_is_last,
                                                       initial_cache_element,
                                                       prg_info);
    return cache_element;
}


KmerIndexCache initial_kmer_index_cache(const Pattern &full_kmer,
                                        const PRG_Info &prg_info) {
    KmerIndexCache cache;

    for (auto it = full_kmer.rbegin(); it != full_kmer.rend(); ++it) {
        const auto &base = *it;
        const bool kmer_base_is_last = it == full_kmer.rbegin();

        if (cache.empty()) {
            auto cache_element = get_initial_cache_element(base, prg_info);
            cache.emplace_back(cache_element);
            continue;
        }

        const auto &last_cache_element = cache.back();
        if (last_cache_element.search_states.empty()) {
            cache.emplace_back(CacheElement {});
            continue;
        }

        const auto &new_cache_element = get_next_cache_element(base,
                                                               kmer_base_is_last,
                                                               last_cache_element,
                                                               prg_info);
        cache.emplace_back(new_cache_element);
    }
    return cache;
}


void update_kmer_index_cache(KmerIndexCache &cache,
                             const Pattern &kmer_prefix_diff,
                             const int kmer_size,
                             const PRG_Info &prg_info) {
    if (kmer_prefix_diff.size() == kmer_size) {
        auto &full_kmer = kmer_prefix_diff;
        cache = initial_kmer_index_cache(full_kmer, prg_info);
        return;
    }

    const auto truncated_cache_size = kmer_size - kmer_prefix_diff.size();
    cache.resize(truncated_cache_size);

    for (auto it = kmer_prefix_diff.rbegin(); it != kmer_prefix_diff.rend(); ++it) {
        const auto &base = *it;
        // the last kmer base is only handled by initial_kmer_index_cache(.)
        const bool kmer_base_is_last = false;

        auto &last_cache_element = cache.back();
        const auto &new_cache_element = get_next_cache_element(base,
                                                               kmer_base_is_last,
                                                               last_cache_element,
                                                               prg_info);
        cache.emplace_back(new_cache_element);
    }
}


void update_full_kmer(Pattern &full_kmer,
                      const Pattern &kmer_prefix_diff,
                      const int kmer_size) {
    if (kmer_prefix_diff.size() == kmer_size) {
        full_kmer = kmer_prefix_diff;
        return;
    }

    auto start_idx = 0;
    for (const auto &base: kmer_prefix_diff)
        full_kmer[start_idx++] = base;
}


KmerIndex index_kmers(const Patterns &kmer_prefix_diffs,
                      const int kmer_size,
                      const PRG_Info &prg_info) {
    KmerIndex kmer_index;
    KmerIndexCache cache;
    Pattern full_kmer;

    auto total_num_kmers = kmer_prefix_diffs.size();
    std::cout << "Total number of unique kmers: "
              << total_num_kmers
              << std::endl << std::endl;

    auto count = 0;
    for (const auto &kmer_prefix_diff: kmer_prefix_diffs) {
        if (count > 0 and count % 50000 == 0)
            std::cout << "Progress: "
                      << count << " of " << total_num_kmers
                      << std::endl;
        count++;

        update_full_kmer(full_kmer,
                         kmer_prefix_diff,
                         kmer_size);

        update_kmer_index_cache(cache,
                                kmer_prefix_diff,
                                kmer_size,
                                prg_info);

        const auto &last_cache_element = cache.back();
        if (not last_cache_element.search_states.empty())
            kmer_index[full_kmer] = last_cache_element.search_states;
    }
    return kmer_index;
}


void generate_kmer_index(const Parameters &parameters,
                         const PRG_Info &prg_info) {
    Patterns kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                       prg_info);
    KmerIndex kmer_index = index_kmers(kmer_prefix_diffs, parameters.kmers_size, prg_info);
    dump_kmer_index(kmer_index, parameters);
}


static inline std::string &left_trim(std::string &s) {
    //trim from start
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}


static inline std::string &right_trim(std::string &s) {
    // trim from end
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}


static inline std::string &trim(std::string &s) {
    // trim from both ends
    return left_trim(right_trim(s));
}


std::vector<std::string> split(const std::string &cad,
                               const std::string &delim) {
    std::vector<std::string> v;
    uint64_t p, q, d;

    q = cad.size();
    p = 0;
    d = strlen(delim.c_str());
    while (p < q) {
        uint64_t posfound = cad.find(delim, p);
        std::string token;

        if (posfound >= 0)
            token = cad.substr(p, posfound - p);
        else
            token = cad.substr(p, cad.size() + 1);
        p += token.size() + d;
        trim(token);
        v.push_back(token);
    }
    return v;
}


std::vector<SA_Interval> parse_sa_intervals(const std::string &full_sa_intervals_str) {
    std::vector<std::string> split_sa_intervals = split(full_sa_intervals_str, " ");
    std::vector<SA_Interval> sa_intervals;
    for (auto i = 0; i < split_sa_intervals.size(); i = i + 2) {
        auto sa_interval_start = (uint64_t) std::stoi(split_sa_intervals[i]);
        auto sa_interval_end = (uint64_t) std::stoi(split_sa_intervals[i + 1]);
        const SA_Interval &sa_interval = std::make_pair(sa_interval_start, sa_interval_end);
        sa_intervals.emplace_back(sa_interval);
    }
    return sa_intervals;
}


VariantSitePath parse_variant_site_path(const std::string &path_data_str) {
    const auto path_split = split(path_data_str, " ");
    if (path_split.empty())
        return VariantSitePath {};

    VariantSitePath variant_site_path;
    for (auto i = 0; i < path_split.size() - 1; i += 2) {
        auto variant_site_marker = (Marker) std::stoi(path_split[i]);
        auto allele_id = (AlleleId) std::stoi(path_split[i + 1]);
        VariantSite variant_site = {variant_site_marker, allele_id};
        variant_site_path.emplace_back(variant_site);
    }
    return variant_site_path;
}


VariantSitePaths parse_vairant_site_paths(const std::vector<std::string> &index_entry_parts) {
    VariantSitePaths variant_site_paths;
    for (uint64_t i = 1; i < index_entry_parts.size(); i++) {
        const auto &path_data_str = index_entry_parts[i];
        const auto &variant_site_path = parse_variant_site_path(path_data_str);
        variant_site_paths.emplace_back(variant_site_path);
    }
    return variant_site_paths;
}


SearchStates parse_kmer_index_entry(const std::string &line) {
    const std::vector<std::string> &parts = split(line, "|");
    const auto &sa_intervals = parse_sa_intervals(parts[0]);
    const auto &variant_site_paths = parse_vairant_site_paths(parts);

    auto i = 0;
    SearchStates search_states;
    for (const auto &variant_site_path: variant_site_paths) {
        const auto &sa_interval = sa_intervals[i++];
        const SearchState search_state = {
                sa_interval,
                variant_site_path
        };
        search_states.emplace_back(search_state);
    }
    return search_states;
}


Pattern deserialize_next_kmer(const uint64_t &kmer_start_index,
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


IndexedKmerStats deserialize_next_stats(const uint64_t &stats_index,
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
        search_states.emplace_back(SearchState {});
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


void parse_sa_intervals(KmerIndex &kmer_index,
                        const sdsl::int_vector<3> &all_kmers,
                        const sdsl::int_vector<> &kmers_stats,
                        const Parameters &parameters) {
    uint64_t sa_interval_index = 0;
    sdsl::int_vector<> sa_intervals;
    sdsl::load_from_file(sa_intervals, parameters.sa_intervals_fpath);

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

            VariantSite site = {marker, allele_id};
            search_state.variant_site_path.emplace_back(site);
        }
    }
}


void parse_paths(KmerIndex &kmer_index,
                 const sdsl::int_vector<3> &all_kmers,
                 const sdsl::int_vector<> &kmers_stats,
                 const Parameters &parameters) {
    uint64_t paths_index = 0;
    sdsl::int_vector<> paths;
    sdsl::load_from_file(paths, parameters.paths_fpath);

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


KmerIndex parse_kmer_index(const Parameters &parameters) {
    KmerIndex kmer_index;

    sdsl::int_vector<3> all_kmers;
    sdsl::load_from_file(all_kmers, parameters.kmers_fpath);

    sdsl::int_vector<> kmers_stats;
    sdsl::load_from_file(kmers_stats, parameters.kmers_stats_fpath);

    parse_sa_intervals(kmer_index, all_kmers, kmers_stats, parameters);
    parse_paths(kmer_index, all_kmers, kmers_stats, parameters);
    return kmer_index;
}


KmerIndex load_kmer_index(const Parameters &parameters) {
    KmerIndex kmer_index;

    uint64_t kmer_start_index = 0;
    sdsl::int_vector<3> all_kmers;
    sdsl::load_from_file(all_kmers, parameters.kmers_fpath);

    std::ifstream fhandle;
    fhandle.open(parameters.kmer_index_fpath);
    std::string line;
    while (std::getline(fhandle, line)) {
        auto kmer = deserialize_next_kmer(kmer_start_index,
                                          all_kmers,
                                          parameters.kmers_size);
        kmer_start_index += parameters.kmers_size;

        auto search_states = parse_kmer_index_entry(line);
        kmer_index[kmer] = search_states;
    }
    return kmer_index;
}


std::string serialize_kmer_index_cache_element(const CacheElement &cache_element) {
    std::stringstream ss;
    ss << "****** Cache Element ******" << std::endl;

    ss << "Base: " << (int) cache_element.base << std::endl;
    ss << "Number of search states: "
       << cache_element.search_states.size()
       << std::endl;

    /*
    for (const auto &search_state: cache_element.search_states) {
        ss << search_state;
    }
     */

    ss << "****** END Cache Element ******" << std::endl;
    return ss.str();
}


std::ostream &operator<<(std::ostream &os, const CacheElement &cache_element) {
    os << serialize_kmer_index_cache_element(cache_element);
    return os;
}


std::string serialize_kmer_index_cache(const KmerIndexCache &cache) {
    std::stringstream ss;
    ss << "****** Cache ******" << std::endl;
    for (const auto &cache_element: cache) {
        ss << cache_element << std::endl;
    }
    ss << "****** END Cache ******" << std::endl;
    return ss.str();
}


std::ostream &operator<<(std::ostream &os, const KmerIndexCache &cache) {
    os << serialize_kmer_index_cache(cache);
    return os;
}
