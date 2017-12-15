#include "parameters.hpp"


#ifndef GRAMTOOLS_COVERAGE_ANALYSIS_HPP
#define GRAMTOOLS_COVERAGE_ANALYSIS_HPP

using AlleleCoverage = std::vector<std::vector<uint32_t>>;

struct QuasimapReadsStats {
    uint64_t all_reads_count = 0;
    uint64_t skipped_reads_count = 0;
    uint64_t mapped_reads_count = 0;
};

QuasimapReadsStats quasimap_reads(const Parameters &parameters,
                             const KmerIndex &kmer_index,
                             const PRG_Info &prg_info);

void quasimap_forward_reverse(QuasimapReadsStats &quasimap_reads_stats,
                              AlleleCoverage &allele_coverage,
                              const Pattern &read,
                              const Parameters &parameters,
                              const KmerIndex &kmer_index,
                              const PRG_Info &prg_info);

bool quasimap_read(const Pattern &read,
                   AlleleCoverage &allele_coverage,
                   const KmerIndex &kmer_index,
                   const PRG_Info &prg_info,
                   const Parameters &parameters);

void record_read_coverage(AlleleCoverage &allele_coverage,
                          const SearchStates &search_states);

void dump_allele_coverage(const AlleleCoverage &allele_coverage,
                          const Parameters &parameters);

AlleleCoverage generate_allele_coverage_structure(const PRG_Info &prg_info);

Pattern get_kmer_from_read(const uint32_t kmer_size, const Pattern &read);

#endif //GRAMTOOLS_COVERAGE_ANALYSIS_HPP
