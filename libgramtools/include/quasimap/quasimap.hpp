#include "parameters.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_QUASIMAP_HPP
#define GRAMTOOLS_QUASIMAP_HPP

struct QuasimapReadsStats {
    uint64_t all_reads_count = 0;
    uint64_t skipped_reads_count = 0;
    uint64_t mapped_reads_count = 0;
};

QuasimapReadsStats quasimap_reads(const Parameters &parameters,
                                  const KmerIndex &kmer_index,
                                  const PRG_Info &prg_info);

void quasimap_forward_reverse(QuasimapReadsStats &quasimap_reads_stats,
                              Coverage &coverage,
                              const Pattern &read,
                              const Parameters &parameters,
                              const KmerIndex &kmer_index,
                              const PRG_Info &prg_info);

bool quasimap_read(const Pattern &read, Coverage &coverage, const KmerIndex &kmer_index, const PRG_Info &prg_info,
                   const Parameters &parameters, const uint32_t &random_seed = 0);

Pattern get_kmer_from_read(const uint32_t& kmer_size, const Pattern &read);

#endif //GRAMTOOLS_QUASIMAP_HPP
