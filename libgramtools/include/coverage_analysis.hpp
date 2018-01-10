#include "parameters.hpp"


#ifndef GRAMTOOLS_COVERAGE_ANALYSIS_HPP
#define GRAMTOOLS_COVERAGE_ANALYSIS_HPP


using AlleleSumCoverage = std::vector<std::vector<uint64_t>>;

template<typename SEQUENCE, typename T>
using sequence_map = std::unordered_map<SEQUENCE, T, seq_hash<SEQUENCE>>;
using AlleleIds = std::vector<AlleleId>;
using GroupedAlleleCounts = sequence_map<AlleleIds, uint64_t>;
using SitesGroupedAlleleCounts = std::vector<GroupedAlleleCounts>;

using BaseCoverage = std::vector<uint64_t>;
using AlleleCoverage = std::vector<BaseCoverage>;
using SitesAlleleBaseCoverage = std::vector<AlleleCoverage>;

struct Coverage {
    AlleleSumCoverage allele_sum_coverage;
    SitesGroupedAlleleCounts grouped_allele_counts;
    SitesAlleleBaseCoverage allele_base_coverage;
};

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

bool quasimap_read(const Pattern &read,
                   Coverage &coverage,
                   const KmerIndex &kmer_index,
                   const PRG_Info &prg_info,
                   const Parameters &parameters);

void record_read_coverage(Coverage &coverage,
                          const SearchStates &search_states);

void dump_coverage(const Coverage &coverage,
                   const Parameters &parameters);

Coverage generate_coverage_structure(const PRG_Info &prg_info);

AlleleSumCoverage generate_allele_sum_coverage_structure(const PRG_Info &prg_info);

SitesAlleleBaseCoverage generate_base_coverage_structure(const PRG_Info &prg_info);

Pattern get_kmer_from_read(const uint32_t kmer_size, const Pattern &read);

#endif //GRAMTOOLS_COVERAGE_ANALYSIS_HPP
