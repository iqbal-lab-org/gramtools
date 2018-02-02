#include "common/utils.hpp"


#ifndef GRAMTOOLS_COVERAGE_TYPES_HPP
#define GRAMTOOLS_COVERAGE_TYPES_HPP

using AlleleSumCoverage = std::vector<std::vector<uint64_t>>;

template<typename SEQUENCE, typename T>
using SequenceHashMap = std::unordered_map<SEQUENCE, T, sequence_hash<SEQUENCE>>;
using AlleleIds = std::vector<AlleleId>;
using GroupedAlleleCounts = SequenceHashMap<AlleleIds, uint64_t>;
using SitesGroupedAlleleCounts = std::vector<GroupedAlleleCounts>;
using AlleleGroupHash = SequenceHashMap<AlleleIds, uint64_t>;

using BaseCoverage = std::vector<uint64_t>;
using AlleleCoverage = std::vector<BaseCoverage>;
using SitesAlleleBaseCoverage = std::vector<AlleleCoverage>;

struct Coverage {
    AlleleSumCoverage allele_sum_coverage;
    SitesGroupedAlleleCounts grouped_allele_counts;
    SitesAlleleBaseCoverage allele_base_coverage;
};

#endif //GRAMTOOLS_COVERAGE_TYPES_HPP
