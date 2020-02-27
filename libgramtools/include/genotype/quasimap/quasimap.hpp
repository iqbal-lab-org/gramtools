/** @file
 * Loads the data structures supporting vBWT search, and maps reads to the prg.
 */

#ifndef GRAMTOOLS_QUASIMAP_HPP
#define GRAMTOOLS_QUASIMAP_HPP

#include "genotype/parameters.hpp"
#include "sequence_read/seqread.hpp"

#include "common/utils.hpp"
#include "genotype/read_stats.hpp"
#include "build/kmer_index/kmer_index_types.hpp"
#include "coverage/common.hpp"

#include "search/encapsulated_search.hpp"

namespace gram {


    struct QuasimapReadsStats {
        uint64_t all_reads_count = 0;
        uint64_t skipped_reads_count = 0;
        uint64_t mapped_reads_count = 0;
        Coverage coverage = {};
    };

    /**
     * For each read file, quasimap reads.
     */
    QuasimapReadsStats quasimap_reads(const GenotypeParams &parameters,
                                      const KmerIndex &kmer_index,
                                      const PRG_Info &prg_info,
                                      ReadStats &readstats);

    /**
     * Load and process (ie map) reads from a given read file using a buffer to reduce disk I/O calls
     */
    void
    handle_read_file(QuasimapReadsStats &quasimap_stats, const std::string &reads_fpath, const GenotypeParams &parameters,
                     const KmerIndex &kmer_index, const PRG_Info &prg_info);

    /**
     * Calls quasimapping routine on a given read (forward mapping), and its reverse complement (reverse mapping)
     */
    void
    quasimap_forward_reverse(QuasimapReadsStats &quasimap_stats, const Sequence &read, const GenotypeParams &parameters,
                             const KmerIndex &kmer_index, const PRG_Info &prg_info);

    /**
     * Map a read to the prg, starting from the precomputed set of search states using the rightmost kmer in the read.
     * @param coverage object in which mapping statistics are recorded.
     * @param kmer_index object holding the pre-computed mappings for kmers. the first kmer in the read will be seeded this way.
     * @param prg_info object holding all data structures necessary for vBWT, including `gram::FM_Index`.
     * @return
     */
    bool quasimap_read(const Sequence &read, Coverage &coverage, const KmerIndex &kmer_index, const PRG_Info &prg_info,
                       const GenotypeParams &parameters);

    Sequence get_kmer_from_read(const uint32_t &kmer_size, const Sequence &read);

    /**
     * Generates a list of `SearchState`s from a read and a kmer, which is 3'-most kmer in the read.
     * The kmer_index is queried to generate an initial set of `SearchState`s (precomputed at `build` stage) to start from.
     * @return SearchStates: a list of `SearchState`s, which at core are an SA interval and a path through the prg (marker-allele ID pairs)
     */
    SearchStates search_read_backwards(const Sequence &read,
                                       const Sequence &kmer,
                                       const KmerIndex &kmer_index,
                                       const PRG_Info &prg_info);

    /**
    * **The key read mapping procedure**.
    * First updates SA_intervals to search next based on variant marker presence.
    * Then executes regular backward search.
    */
    SearchStates process_read_char_search_states(const int_Base &pattern_char,
                                                 const SearchStates &old_search_states,
                                                 const PRG_Info &prg_info);



    Sequence reverse_complement_read(const Sequence &read);
}
#endif //GRAMTOOLS_QUASIMAP_HPP
