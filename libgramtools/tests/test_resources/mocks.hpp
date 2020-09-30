#include "common/random.hpp"
#include "genotype/infer/types.hpp"
#include "genotype/read_stats.hpp"
#include "gmock/gmock.h"

using namespace gram;

class MockRandomGenerator : public gram::RandomGenerator {
 public:
  MOCK_METHOD(uint32_t, generate, (uint32_t, uint32_t), (override));
};

class MockReadStats : public gram::AbstractReadStats {
 public:
  MOCK_METHOD(allele_and_cov, extract_max_coverage_allele,
              (SitesGroupedAlleleCounts const&, covG_ptr, covG_ptr),
              (override));
};
