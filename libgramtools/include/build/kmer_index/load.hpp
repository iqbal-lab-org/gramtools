/** @file
 * Routines for restoring (deserialising) a `gram::KmerIndex` from a set of files describing it.
 * @see dump.hpp for an explanation of the serialising procedure.
 */
#include "common/utils.hpp"
#include "common/parameters.hpp"

#include "build.hpp"
#include "kmer_index_types.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_LOAD_HPP
#define GRAMTOOLS_KMER_INDEX_LOAD_HPP

namespace gram {

    /**
     * Extracts a nucleotide kmer from an integer vector representing all indexed kmers.
     */
    Sequence deserialize_next_kmer(const uint64_t &kmer_start_index,
                                   const sdsl::int_vector<3> &all_kmers,
                                   const uint32_t &kmers_size);

    /**
     * Extracts the statistics of a single indexed kmer.
     * @return `gram::IndexedKmerStats` giving the number of `gram::SearchState`s and the length of each of their associated `gram::variant_site_path`.
     */
    IndexedKmerStats deserialize_next_stats(const uint64_t &stats_index,
                                            const sdsl::int_vector<> &kmers_stats);

    /**
     * Rebuilds `gram::SearchStates` for each indexed kmer, populating each `gram::SearchState` with an `gram::SA_Interval`.
     */
    void parse_sa_intervals(KmerIndex &kmer_index,
                            const sdsl::int_vector<3> &all_kmers,
                            const sdsl::int_vector<> &kmers_stats,
                            CommonParameters const &parameters);

    void parse_paths(KmerIndex &kmer_index,
                     const sdsl::int_vector<3> &all_kmers,
                     const sdsl::int_vector<> &kmers_stats,
                     CommonParameters const &parameters);

    namespace kmer_index {
        /**
         * Rebuild a `gram::KmerIndex` from serialised file in a gramtools `build` produced directory.
         */
        KmerIndex load(CommonParameters const &parameters);
    }

}

#endif //GRAMTOOLS_KMER_INDEX_LOAD_HPP
