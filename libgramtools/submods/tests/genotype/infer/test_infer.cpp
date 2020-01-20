#include "gtest/gtest.h"
#include "genotype/infer/genotyping_models.hpp"
#include "genotype/infer/infer.hpp"

class TestLevelGenotyperModel : public ::testing::Test {
protected:
    allele_vector alleles{
        Allele{ "ATCACC", {0, 0, 1, 1, 2, 2}, 0 },
        Allele{ "GGGCC", {10, 12, 12, 14, 14}, 1 },
    };

    GroupedAlleleCounts gp_counts{
            { {0}, 1 },
            { {0, 1}, 1 }, // The allele sequences support one read being able to map like this
            { {1}, 13 },
    };

    double mean_cov_depth{15}, mean_pb_error{0.01};
    likelihood_related_stats l_stats = LevelGenotyper::make_l_stats(mean_cov_depth, mean_pb_error);
};

TEST_F(TestLevelGenotyperModel, Given0MeanCoverage_ReturnsNullGenotypedSite) {
    l_stats.mean_cov_depth = 0;
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats);

    EXPECT_TRUE(genotyped.get_site()->is_null());
}
