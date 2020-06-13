/** @file
 * Defines coverage related operations for allele sum coverage.
 * `AlleleSumCoverage` stores the sum of all reads mapped for each allele of
 * each variant site.
 */
#include "coverage_common.hpp"
#include "genotype/quasimap/coverage/types.hpp"
#include "genotype/quasimap/search/types.hpp"

#ifndef GRAMTOOLS_ALLELE_SUM_HPP
#define GRAMTOOLS_ALLELE_SUM_HPP

namespace gram::coverage {
namespace generate {
/**
 * Generates the coverage structure recording allele sum counts.
 * Iterates over the bubbles of the coverage Graph to do so.
 * @return `gram::AlleleSumCoverage` A vector of vectors of integers. The
 * top-level vector represents each site, as a vector of allele counts.
 */
AlleleSumCoverage allele_sum_structure(const PRG_Info &prg_info);
}  // namespace generate

namespace record {
/**
 * Increments each site/allele combination compatible with the mapped read.
 * @param coverage The `Coverage` structure common to all mapped reads.
 * @param compatible_loci The selected `SearchStates` for recording coverage.
 */
void allele_sum(Coverage &coverage, const uniqueLoci &compatible_loci);
}  // namespace record

namespace dump {
void allele_sum(const Coverage &coverage, const GenotypeParams &parameters);
}
}  // namespace gram::coverage

#endif  // GRAMTOOLS_ALLELE_SUM_HPP
