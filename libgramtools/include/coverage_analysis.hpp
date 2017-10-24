#include "parameters.hpp"


#ifndef GRAMTOOLS_COVERAGE_ANALYSIS_HPP
#define GRAMTOOLS_COVERAGE_ANALYSIS_HPP

using AlleleCoverage = std::vector<std::vector<uint32_t>>;

uint64_t quasimap_reads(const Parameters &params,
                        const KmerIndex &kmer_index,
                        const PRG_Info &prg_info);

bool quasimap_read(const Pattern &read,
                   AlleleCoverage &allele_coverage,
                   const KmerIndex &kmer_index,
                   const PRG_Info &prg_info,
                   const Parameters &params);

void record_read_coverage(AlleleCoverage &allele_coverage,
                          const SearchStates &search_states);

void dump_allele_coverage(const AlleleCoverage &allele_coverage,
                          const Parameters &params);

AlleleCoverage make_allele_coverage_structure(const PRG_Info &prg_info);

Pattern get_kmer_from_read(const uint32_t kmer_size, const Pattern &read);

#endif //GRAMTOOLS_COVERAGE_ANALYSIS_HPP
