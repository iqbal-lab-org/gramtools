/** @file
 * Defines coverage related types.
 */

#ifndef GRAMTOOLS_COVERAGE_TYPES_HPP
#define GRAMTOOLS_COVERAGE_TYPES_HPP

#include "prg/types.hpp"
#include "common/utils.hpp"
#include "genotype/parameters.hpp"

namespace gram {
    using AlleleSumCoverage = std::vector<PerAlleleCoverage>; /**<Number of reads mapped per allele for each variant site.*/

    /** Vector of `gram::AlleleId`. Used to store different alleles of the same variant site both by a read.*/
    using AlleleIds = std::vector<AlleleId>;
    /* An unordered_map associating a group of alleles (`gram::AlleleIds`) with a count of how many reads mapped to this group.*/
    using GroupedAlleleCounts = SequenceHashMap<AlleleIds, CovCount>;
    /** A vector containing unordered_maps of allele group counts.
     * There is one such map per variant site in the prg.*/
    using SitesGroupedAlleleCounts = std::vector<GroupedAlleleCounts>;

    using AlleleGroupHash = SequenceHashMap<AlleleIds, uint64_t>;

    using SitePbCoverage = std::vector<PerBaseCoverage>; /**< `gram::PerBaseCoverage` for each allele of a variant site. */
    using SitesAlleleBaseCoverage = std::vector<SitePbCoverage>; /**< Vector of gram::AlleleCoverage, one for each variant site in the prg. */

    /**
     * Groups together all coverage metrics to record.
     */
    struct Coverage {
        AlleleSumCoverage allele_sum_coverage;
        SitesGroupedAlleleCounts grouped_allele_counts;
        SitesAlleleBaseCoverage allele_base_coverage;
    };
}

#endif //GRAMTOOLS_COVERAGE_TYPES_HPP
