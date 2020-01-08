#include <algorithm>
#include <thread>
#include <unordered_map>

#include "genotype/quasimap/search/vBWT_jump.hpp"
#include "genotype/quasimap/search/BWT_search.hpp"
#include "kmer_index/load.hpp"
#include "kmer_index/kmers.hpp"
#include "kmer_index/build.hpp"

using namespace gram;


/**
 * Variant aware backward search on the next base of a kmer.
 * @see process_markers_search_states()
 * @see search_base_backwards()
 */
CacheElement get_next_cache_element(const int_Base &base,
                                    const bool kmer_base_is_first_processed,
                                    const CacheElement &last_cache_element,
                                    const PRG_Info &prg_info) {
    const auto &old_search_states = last_cache_element.search_states;
    SearchStates new_search_states;
    if (not kmer_base_is_first_processed) {
        new_search_states = process_markers_search_states(old_search_states,
                                                          prg_info);
    } else {
        new_search_states = old_search_states;
    }
    new_search_states = search_base_backwards(base,
                                              new_search_states,
                                              prg_info);
    return CacheElement{
            new_search_states,
            base
    };
}

/**
 * Wrapper around `get_next_cache_element` to deal with first search of a base in the prg.
 * @see get_next_cache_element()
 */
CacheElement get_initial_cache_element(const int_Base &base,
                                       const PRG_Info &prg_info) {
    // Start with the full SA interval, over the whole PRG
    SearchState search_state = {
            SA_Interval{0, prg_info.fm_index.size() - 1}
    };
    SearchStates search_states = {search_state};
    CacheElement full_sa_interval = {search_states};

    bool kmer_base_is_first_processed = true;
    const auto &cache_element = get_next_cache_element(base,
                                                       kmer_base_is_first_processed,
                                                       full_sa_interval,
                                                       prg_info);
    return cache_element;
}


/**
 * Routine for updating `SearchStates` using backward search starting from the `cache`
 * @param cache a `KmerIndexCache`: list of `CacheElement`s, which contain one set of `SearchStates` and one `Base`
 * @param full_kmer: a `Pattern`, which is a vector of `Base`s (`uint8`s)
 */
void build_kmer_cache(KmerIndexCache &cache,
                      const Pattern &kmer_prefix_diff,
                      const int kmer_size,
                      const PRG_Info &prg_info) {
    auto it = kmer_prefix_diff.rbegin();

    if (kmer_prefix_diff.size() == kmer_size) {
    // Case: a fully new kmer is encountered. No reuse of `cache`d search possible. Call a full search on the kmer.
        const auto &base = *it;
        cache.resize(0); // Empty the `cache`
        auto cache_element = get_initial_cache_element(base, prg_info);
        cache.emplace_back(cache_element);
        ++it;
    }

    else {
        const auto preserved_cache_size = kmer_size - kmer_prefix_diff.size();
        cache.resize(preserved_cache_size);
    }

    for (; it != kmer_prefix_diff.rend(); ++it) {
        const auto &base = *it;
        // the right-most kmer base (first processed) is only ever handled by `get_initial_cache_element`
        const bool kmer_base_is_first_processed = false;

        auto &last_cache_element = cache.back();
        const auto &new_cache_element = get_next_cache_element(base,
                                                               kmer_base_is_first_processed,
                                                               last_cache_element,
                                                               prg_info);
        cache.emplace_back(new_cache_element);
    }
}

/**
 * Modify the `full_kmer` using the set of differences recorded previously.
 */
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

KmerIndex gram::index_kmers(const Patterns &kmer_prefix_diffs,
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

        // Obtain the full kmer from the previous kmer and the current prefix_diff
        update_full_kmer(full_kmer,
                         kmer_prefix_diff,
                         kmer_size);

        // Call cache update routine
        build_kmer_cache(cache,
                         kmer_prefix_diff,
                         kmer_size,
                         prg_info);

        // Associate the current kmer with the resulting `SearchStates`, if they are not empty.
        const auto &last_cache_element = cache.back();
        if (not last_cache_element.search_states.empty())
            kmer_index[full_kmer] = last_cache_element.search_states;
    }
    return kmer_index;
}

/**
 * Highest level indexing routine.
 * @see get_kmer_prefix_diffs()
 * @see index_kmers()
 */
KmerIndex gram::kmer_index::build(const Parameters &parameters,
                                  const PRG_Info &prg_info) {
    // Extract all relevant kmers and generate the minimal differences between them.
    Patterns kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                       prg_info);
    std::cout << "Indexing kmers" << std::endl;
    KmerIndex kmer_index = index_kmers(kmer_prefix_diffs, parameters.kmers_size, prg_info);
    return kmer_index;
}
