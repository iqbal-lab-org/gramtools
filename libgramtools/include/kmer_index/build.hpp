#include "common/utils.hpp"
#include "common/parameters.hpp"
#include "prg/prg.hpp"
#include "prg/fm_index.hpp"
#include "search/search_types.hpp"
#include "kmer_index_types.hpp"
#include "kmers.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_BUILD_HPP
#define GRAMTOOLS_KMER_INDEX_BUILD_HPP

struct KmerIndexStats {
    uint64_t count_kmers;
    uint64_t count_search_states;
    uint64_t count_total_path_elements;
};

struct IndexedKmerStats {
    uint64_t count_search_states;
    std::vector<uint64_t> path_lengths;

    bool operator==(const IndexedKmerStats &other) const {
        return this->count_search_states == other.count_search_states
               and this->path_lengths == other.path_lengths;
    };
};

KmerIndex index_kmers(const Patterns &kmers, const int kmer_size, const PRG_Info &prg_info);

namespace kmer_index {
    KmerIndex build(const Parameters &parameters,
                    const PRG_Info &prg_info);
}

#endif //GRAMTOOLS_KMER_INDEX_BUILD_HPP
