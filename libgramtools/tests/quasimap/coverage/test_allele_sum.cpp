#include <cctype>

#include "gtest/gtest.h"

#include "../../test_utils.hpp"
#include "quasimap/coverage/allele_sum.hpp"


using namespace gram;


TEST(AlleleSumCoverage, GivenOneVariantSite_CorrectAlleleSumCoverageStructure) {
    auto prg_raw = "gcgct5gg6agtg5ctgt";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_sum_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSumCoverage, GivenTwoVariantSite_CorrectAlleleSumCoverageStructure) {
    auto prg_raw = "gcgct5gg6agtg5cccc7t8g7t";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_sum_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSumCoverage, GivenThreeVariantSites_CorrectAlleleSumCoverageStructure) {
    auto prg_raw = "5gg6agtg5c7t8g8c7t9ccccc10t9";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_sum_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0},
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}