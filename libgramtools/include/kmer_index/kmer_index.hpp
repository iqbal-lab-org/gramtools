#include "common/utils.hpp"
#include "common/parameters.hpp"
#include "prg/prg.hpp"
#include "prg/fm_index.hpp"
#include "search/search_types.hpp"
#include "kmer_index_types.hpp"
#include "kmers.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_HPP
#define GRAMTOOLS_KMER_INDEX_HPP

struct KmerIndexStats {
    uint64_t count_kmers;
    uint64_t count_search_states;
    uint64_t count_total_path_elements;
};

KmerIndexStats calculate_stats(const KmerIndex &kmer_index);

void dump_kmers_stats(const KmerIndexStats &stats,
                      const sdsl::int_vector<3> &all_kmers,
                      const KmerIndex &kmer_index,
                      const Parameters &parameters);

sdsl::int_vector<3> dump_kmers(const KmerIndex &kmer_index,
                               const Parameters &parameters);

void dump_sa_intervals(const KmerIndexStats &stats,
                       const sdsl::int_vector<3> &all_kmers,
                       const KmerIndex &kmer_index,
                       const Parameters &parameters);

void dump_paths(const KmerIndexStats &stats,
                const sdsl::int_vector<3> &all_kmers,
                const KmerIndex &kmer_index,
                const Parameters &parameters);

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
    void dump(const KmerIndex &kmer_index, const Parameters &parameters);

    KmerIndex build(const Parameters &parameters,
                    const PRG_Info &prg_info);
}

#endif //GRAMTOOLS_KMER_INDEX_HPP
