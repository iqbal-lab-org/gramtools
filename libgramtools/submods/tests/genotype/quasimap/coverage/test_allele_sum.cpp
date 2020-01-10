#include <cctype>

#include "gtest/gtest.h"

#include "src_common/generate_prg.hpp"
#include "genotype/quasimap/coverage/allele_sum.hpp"


using namespace gram;


TEST(AlleleSumCoverage, GivenOneVariantSite_CorrectAlleleSumCoverageStructure) {
    auto prg_raw = encode_prg("gcgct5gg6agtg6ctgt");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_sum_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSumCoverage, GivenTwoVariantSite_CorrectAlleleSumCoverageStructure) {
    auto prg_raw = encode_prg("gcgct5gg6agtg6cccc7t8g8t");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_sum_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSumCoverage, GivenThreeVariantSites_CorrectAlleleSumCoverageStructure) {
    auto prg_raw = encode_prg("5gg6agtg6c7t8g8c8t9ccccc10t10");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_sum_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0},
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSumCoverage_Nested, GivenFourNestedVariantSites_CorrectAlleleSumCoverageStructure) {
    auto prg_raw = prg_string_to_ints("A[GC,T[C,G[CG,C]]A,C]AAT[A,C,G,T]C");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_sum_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 0},
            {0, 0},
            {0, 0, 0, 0}
    };
    EXPECT_EQ(result, expected);
}
