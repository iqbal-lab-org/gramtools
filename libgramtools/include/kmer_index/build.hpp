/**
 * @file
 * Build an index into search states in the PRG in order to speed up read quasimapping.
 * The index maps a kmer to a set of `SearchStates`: variant site paths and SA intervals.
 * Then during quasimapping the search is initialised to the read's last kmer index entry.
 *
 * Kmer size is passed as a parameter. Either all kmers of this size are enumerated, or all kmers associated with variant
 * sites in the PRG are searched for in a given PRG.
 */
#include "common/utils.hpp"
#include "common/parameters.hpp"
#include "prg/prg.hpp"
#include "prg/make_data_structures.hpp"
#include "search/search_types.hpp"
#include "kmer_index_types.hpp"
#include "kmers.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_BUILD_HPP
#define GRAMTOOLS_KMER_INDEX_BUILD_HPP

namespace gram {

    /**
     * Stores the total number of indexed kmers, the total number of computed `gram::SA_Interval`s, and the total number of traversed `gram::VariantLocus`.
     */
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


    /**
     * For each kmer, find its `SearchStates` and populate the `KmerIndex`.
     * @see update_full_kmer()
     * @see update_kmer_index_cache()
     */
    KmerIndex index_kmers(const Patterns &kmers, const int kmer_size, const PRG_Info &prg_info);

    namespace kmer_index {
        KmerIndex build(const Parameters &parameters,
                        const PRG_Info &prg_info);
    }

}

#endif //GRAMTOOLS_KMER_INDEX_BUILD_HPP
