/** @file
 * Defines coverage related types.
 */
#include "common/utils.hpp"


#ifndef GRAMTOOLS_COVERAGE_TYPES_HPP
#define GRAMTOOLS_COVERAGE_TYPES_HPP

namespace gram {
    using AlleleSumCoverage = std::vector<std::vector<uint64_t>>; /**<Number of reads mapped per allele for each variant site.*/

    // TODO: maybe remove this line? As it is used in utils.hpp
    template<typename SEQUENCE, typename T>
    using SequenceHashMap = std::unordered_map<SEQUENCE, T, sequence_hash < SEQUENCE>>;

    /** Vector of `gram::AlleleId`. Used to store different alleles of the same variant site both by a read.*/
    using AlleleIds = std::vector<AlleleId>;
    /* An unordered_map associating a group of alleles (`gram::AlleleIds`) with a count of how many reads mapped to this group.*/
    using GroupedAlleleCounts = SequenceHashMap<AlleleIds, uint64_t>;
    /** A vector containing unordered_maps of allele group counts.
     * There is one such map per variant site in the prg.*/
    using SitesGroupedAlleleCounts = std::vector<GroupedAlleleCounts>;

    using AlleleGroupHash = SequenceHashMap<AlleleIds, uint64_t>;

    using BaseCoverage = std::vector<uint16_t>; /**< Number of reads mapped to each base of an allele */
    using AlleleCoverage = std::vector<BaseCoverage>; /**< `gram::BaseCoverage` for each allele of a variant site. */
    using SitesAlleleBaseCoverage = std::vector<AlleleCoverage>; /**< Vector of gram::AlleleCoverage, one for each variant site in the prg. */

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
