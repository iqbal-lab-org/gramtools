#include <gtest/gtest.h>

#include "genotype/infer/types.hpp"

using namespace gram;

TEST(Alleles, CombineAlleles) {
  Allele a1{"ATA", {0, 1, 0}, 0};
  Allele a2{"TT", {2, 0}, 1};
  auto result = a1 + a2;
  Allele expected{"ATATT", {0, 1, 0, 2, 0}, 0};
  EXPECT_EQ(result, expected);
}

TEST(Alleles, GetAverageCoverage) {
  Allele a1{"ATAT", {2, 5, 0, 3}, 0};
  auto result = a1.get_average_cov();
  EXPECT_DOUBLE_EQ(result, 2.5);
}