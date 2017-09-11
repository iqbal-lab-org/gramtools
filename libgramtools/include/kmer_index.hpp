#include <sdsl/suffix_trees.hpp>
#include "fm_index.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_HPP
#define GRAMTOOLS_KMER_INDEX_HPP

using Kmer = std::vector<uint8_t>;
using Kmers = std::vector<Kmer>;

using SA_Interval = std::pair<uint64_t, uint64_t>;
using SA_Intervals = std::list<SA_Interval>;

struct CacheElement {
    SA_Intervals sa_intervals;
    Sites sites;
    uint8_t base = 0;
};

using KmerIndexCache = std::list<CacheElement>;

void update_kmer_index_cache(KmerIndexCache &cache,
                             const Kmer &kmer_suffix_diff,
                             const uint64_t kmer_size,
                             const PRG_Info &prg_info);

void generate_kmer_index(const Kmers &kmer_suffix_diffs, const PRG_Info &prg_info);


#endif //GRAMTOOLS_KMER_INDEX_HPP
