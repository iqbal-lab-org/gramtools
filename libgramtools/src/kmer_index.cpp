#include <list>

#include "prg.hpp"
#include "bidir_search_bwd.hpp"
#include "kmer_index.hpp"


void update_kmer_index_cache(KmerIndexCache &cache,
                             const Kmer &kmer_suffix_diff,
                             const uint64_t kmer_size,
                             const PRG_Info &prg_info) {

    // truncate cache
    const auto new_cache_size = kmer_size - kmer_suffix_diff.size();
    cache.resize(new_cache_size);

    for (const auto &base: kmer_suffix_diff) {
        CacheElement new_cache_element;
        new_cache_element.base = base;

        if (cache.size() == 0) {
            new_cache_element.sa_intervals = {{0, prg_info.fm_index.size()}};
            new_cache_element.sites = {Site()};
        } else {
            const auto &last_cache_element = cache.back();
            new_cache_element.sa_intervals = last_cache_element.sa_intervals;
            new_cache_element.sites = last_cache_element.sites;
        }

        reduce_sa_intervals(base,
                            new_cache_element.sa_intervals,
                            new_cache_element.sites,
                            false, true,
                            prg_info.allele_mask,
                            prg_info.max_alphabet_num,
                            false,
                            prg_info.dna_rank,
                            prg_info.fm_index);

        cache.emplace_back(new_cache_element);
    }
}


void generate_kmer_index(const Kmers &kmer_suffix_diffs, const PRG_Info &prg_info) {
    const uint64_t kmer_size = kmer_suffix_diffs[0].size();
    KmerIndexCache cache;

    for (const auto &kmer_suffix_diff: kmer_suffix_diffs) {
        update_kmer_index_cache(cache,
                                kmer_suffix_diff,
                                kmer_size,
                                prg_info);

        // operate on last cache element
    }
}
