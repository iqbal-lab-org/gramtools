#include "kmer_index/build.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_DUMP_HPP
#define GRAMTOOLS_KMER_INDEX_DUMP_HPP

KmerIndexStats calculate_stats(const KmerIndex &kmer_index);

sdsl::int_vector<3> dump_kmers(const KmerIndex &kmer_index,
                               const Parameters &parameters);

void dump_kmers_stats(const KmerIndexStats &stats,
                      const sdsl::int_vector<3> &all_kmers,
                      const KmerIndex &kmer_index,
                      const Parameters &parameters);

void dump_sa_intervals(const KmerIndexStats &stats,
                       const sdsl::int_vector<3> &all_kmers,
                       const KmerIndex &kmer_index,
                       const Parameters &parameters);

void dump_paths(const KmerIndexStats &stats,
                const sdsl::int_vector<3> &all_kmers,
                const KmerIndex &kmer_index,
                const Parameters &parameters);

namespace kmer_index {
    void dump(const KmerIndex &kmer_index, const Parameters &parameters);
}

#endif //GRAMTOOLS_KMER_INDEX_DUMP_HPP
