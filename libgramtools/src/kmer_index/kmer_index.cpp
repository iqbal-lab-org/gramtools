#include <algorithm>
#include <thread>
#include <unordered_map>
#include <kmer_index/kmers.hpp>

#include "kmer_index/kmer_index.hpp"
#include "search/search.hpp"


std::string dump_kmer(const Pattern &kmer) {
    std::stringstream stream;
    for (const auto &base: kmer) {
        stream << (int) base;
        bool last_base = (&base == &kmer.back());
        if (!last_base)
            stream << ' ';
    }
    return stream.str();
}


std::string dump_sa_intervals(const SearchStates &search_states) {
    std::stringstream stream;

    for (const auto &search_state: search_states) {
        const auto &sa_interval = search_state.sa_interval;
        stream << sa_interval.first
               << " "
               << sa_interval.second;
        bool last_search_state = (&search_state == &search_states.back());
        if (not last_search_state)
            stream << " ";
    }
    return stream.str();
}


std::string dump_variant_site_paths(const SearchStates &search_states) {
    std::stringstream stream;

    for (const auto &search_state: search_states) {
        const auto &variant_site_path = search_state.variant_site_path;

        for (const auto &variant_site: variant_site_path) {
            const auto &site_marker = variant_site.first;
            const auto &allele_id = variant_site.second;

            stream << site_marker << " " << allele_id;

            const bool last_variant_site = &variant_site == &variant_site_path.back();
            if (not last_variant_site)
                stream << " ";
        }
        stream << "|";
    }
    return stream.str();
}


std::string dump_kmer_index_entry(const SearchStates &search_states) {
    std::stringstream stream;
    stream << dump_sa_intervals(search_states) << "|";
    stream << dump_variant_site_paths(search_states);
    return stream.str();
}


void dump_kmer_index(const KmerIndex &kmer_index,
                     const Parameters &parameters) {
    std::ofstream kmer_index_file;
    kmer_index_file.open(parameters.kmer_index_fpath);

    uint64_t i = 0;
    sdsl::int_vector<3> all_kmers(kmer_index.size() * parameters.kmers_size);

    for (const auto &kmer_search_states: kmer_index) {
        const auto &kmer = kmer_search_states.first;
        for (const auto &base: kmer)
            all_kmers[i++] = base;

        const auto &search_states = kmer_search_states.second;
        const auto &kmer_entry = dump_kmer_index_entry(search_states);
        kmer_index_file << kmer_entry << std::endl;
    }

    sdsl::store_to_file(all_kmers, parameters.kmers_fpath);
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
    assert(kmer_start_index <= all_kmers.size() - kmers_size);
    Pattern kmer;
    kmer.reserve(kmers_size);
    for (uint64_t i = kmer_start_index; i < kmer_start_index + kmers_size; ++i)
        kmer.emplace_back(all_kmers[i]);
    return kmer;
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
