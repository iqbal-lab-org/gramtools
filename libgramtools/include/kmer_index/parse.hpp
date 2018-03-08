#include "../common/utils.hpp"
#include "../common/parameters.hpp"

#include "kmer_index.hpp"
#include "kmer_index_types.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_PARSE_HPP
#define GRAMTOOLS_KMER_INDEX_PARSE_HPP

Pattern deserialize_next_kmer(const uint64_t &kmer_start_index,
                              const sdsl::int_vector<3> &all_kmers,
                              const uint32_t &kmers_size);

IndexedKmerStats deserialize_next_stats(const uint64_t &stats_index,
                                        const sdsl::int_vector<> &kmers_stats);

void parse_sa_intervals(KmerIndex &kmer_index,
                        const sdsl::int_vector<3> &all_kmers,
                        const sdsl::int_vector<> &kmers_stats,
                        const Parameters &parameters);

void parse_paths(KmerIndex &kmer_index,
                 const sdsl::int_vector<3> &all_kmers,
                 const sdsl::int_vector<> &kmers_stats,
                 const Parameters &parameters);

namespace kmer_index {
    KmerIndex parse(const Parameters &parameters);
}


#endif //GRAMTOOLS_KMER_INDEX_PARSE_HPP
