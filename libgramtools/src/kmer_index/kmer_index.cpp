#include <algorithm>
#include <thread>
#include <unordered_map>

#include "search/search.hpp"
#include "kmer_index/load.hpp"
#include "kmer_index/kmers.hpp"
#include "kmer_index/kmer_index.hpp"


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


KmerIndex kmer_index::build(const Parameters &parameters,
                            const PRG_Info &prg_info) {
    Patterns kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                       prg_info);
    KmerIndex kmer_index = index_kmers(kmer_prefix_diffs, parameters.kmers_size, prg_info);
    return kmer_index;
}
