/** @file
 * Loads the data structures supporting vBWT search, and maps reads to the prg.
 */
#include "parameters.hpp"
#include "kmer_index/kmer_index_types.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_QUASIMAP_HPP
#define GRAMTOOLS_QUASIMAP_HPP

namespace gram {

    namespace commands::quasimap {
        void run(const Parameters &parameters);
    }

    struct QuasimapReadsStats {
        uint64_t all_reads_count = 0;
        uint64_t skipped_reads_count = 0;
        uint64_t mapped_reads_count = 0;
    };

    /**
     * For each read file, quasimap reads.
     */
    QuasimapReadsStats quasimap_reads(const Parameters &parameters,
                                      const KmerIndex &kmer_index,
                                      const PRG_Info &prg_info);

    /**
     * Load and process (ie map) reads from a given read file using a buffer to reduce disk I/O calls
     */
    void handle_read_file(QuasimapReadsStats &quasimap_stats, Coverage &coverage, const std::string &reads_fpath,
                          const Parameters &parameters, const KmerIndex &kmer_index, const PRG_Info &prg_info);

    /**
     * Calls quasimapping routine on a given read (forward mapping), and its reverse complement (reverse mapping)
     */
    void quasimap_forward_reverse(QuasimapReadsStats &quasimap_reads_stats,
                                  Coverage &coverage,
                                  const Pattern &read,
                                  const Parameters &parameters,
                                  const KmerIndex &kmer_index,
                                  const PRG_Info &prg_info);

    /**
     * Map a read to the prg, starting from the precomputed set of search states using the rightmost kmer in the read.
     * @param coverage object in which mapping statistics are recorded.
     * @param kmer_index object holding the pre-computed mappings for kmers. the first kmer in the read will be seeded this way.
     * @param prg_info object holding all data structures necessary for vBWT, including `gram::FM_Index`.
     * @return
     */
    bool quasimap_read(const Pattern &read, Coverage &coverage, const KmerIndex &kmer_index, const PRG_Info &prg_info,
                       const Parameters &parameters);

    Pattern get_kmer_from_read(const uint32_t &kmer_size, const Pattern &read);

}

#endif //GRAMTOOLS_QUASIMAP_HPP
