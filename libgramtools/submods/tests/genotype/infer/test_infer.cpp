#include "gtest/gtest.h"
#include "genotype/infer/genotyping_models.hpp"
#include "genotype/infer/infer.hpp"

/**
 * LevelGenotyperModel tests
 */
TEST(TestLevelGenotyperModel_Failure, GivenOneAlleleOnly_Breaks){
    // No likelihood ratio if only one allele. Note this should not present itself if allele extraction
    // works correctly, as any bubble has at least 2 alleles.
    allele_vector alleles{
        Allele{"ACGT", {1, 1, 1, 1}, 0}
    };
    GroupedAlleleCounts gp_counts;
    likelihood_related_stats l_stats;
    EXPECT_DEATH(LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats), "");
}

class TestLevelGenotyperModel_NullGTs : public ::testing::Test {
protected:
    allele_vector alleles{
            Allele{ "A", {0}, 0 },
            Allele{ "G", {0}, 1 },
    };

    GroupedAlleleCounts gp_counts{};
    double mean_cov_depth{15}, mean_pb_error{0.01};
    likelihood_related_stats l_stats = LevelGenotyper::make_l_stats(mean_cov_depth, mean_pb_error);
};

TEST_F(TestLevelGenotyperModel_NullGTs , Given0MeanCoverage_ReturnsNullGenotypedSite) {
    l_stats.mean_cov_depth = 0;
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats);

    EXPECT_TRUE(genotyped.get_site()->is_null());
}

TEST_F(TestLevelGenotyperModel_NullGTs , GivenNoCoverageOnAllAlleles_ReturnsNullGenotypedSite) {
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats);

    EXPECT_TRUE(genotyped.get_site()->is_null());
}


class TestLevelGenotyperModel_TwoAllelesWithCoverage : public ::testing::Test {
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


TEST_F(TestLevelGenotyperModel_TwoAllelesWithCoverage, GivenCoverage_ReturnsCorrectHaploidCall) {
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats);

    auto gtype = genotyped.get_site()->get_genotype();
    EXPECT_TRUE(std::holds_alternative<GtypedIndices>(gtype));

    GtypedIndices expected_gtype{1};
    EXPECT_EQ(std::get<GtypedIndices>(gtype), expected_gtype);
}


TEST_F(TestLevelGenotyperModel_TwoAllelesWithCoverage, GivenCoverage_ReturnsCorrectDiploidCall) {
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats);

    auto gtype = genotyped.get_site()->get_genotype();
    EXPECT_TRUE(std::holds_alternative<GtypedIndices>(gtype));

    GtypedIndices expected_gtype{1, 1};
    EXPECT_EQ(std::get<GtypedIndices>(gtype), expected_gtype);
}

TEST(TestLevelGenotyperModel_MinosParallel, GivenCoverages_CorrectGenotype){
    // Note : comparing with Minos v0.9.1 commit@7c68ad0 (and with the hom likelihood not halved),
    // which has the same test as this, I get the same likelihoods for the three genotypes
    double mean_cov_depth{20}, mean_pb_error{0.01};

    allele_vector alleles{
            Allele{ "AA", {0, 1}, 0 },
            Allele{ "TT", {20, 19}, 1 },
    };

    GroupedAlleleCounts gp_counts{
            { {0}, 2 },
            { {0, 1}, 1 },
            { {1}, 20 },
    };

    likelihood_related_stats l_stats = LevelGenotyper::make_l_stats(mean_cov_depth, mean_pb_error);

    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats);
    auto gtype = genotyped.get_site()->get_genotype();

    GtypedIndices expected_gtype{1, 1};
    EXPECT_EQ(std::get<GtypedIndices>(gtype), expected_gtype);
}


class TestLevelGenotyperModel_FourAlleles : public ::testing::Test{
   // Simulating a case where each haplogroup has a bubble nested inside it,
   // and those bubbles have been typed heterozygous.
protected:
    allele_vector alleles{
            Allele{ "AATAA", {8, 8, 8, 8, 8}, 0 },
            Allele{ "AAGAA", {7, 7, 7, 7, 7}, 0 },
            Allele{ "GGTGG", {15, 15, 15, 16, 16}, 1 }, // 15 unique + 1 common with next allele
            Allele{ "GGCGG", {14, 14, 14, 15, 15}, 1 }, // 14 unique + 1 common with previous allele
    };

    GroupedAlleleCounts gp_counts{
            { {0}, 15 },
            { {1}, 30 },
    };
    double mean_cov_depth{30}, mean_pb_error{0.01};
    likelihood_related_stats l_stats = LevelGenotyper::make_l_stats(mean_cov_depth, mean_pb_error);
};

TEST_F(TestLevelGenotyperModel_FourAlleles, GivenHaploGroup1SupportingMeanCov_CorrectGenotype){

    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats);
    auto gtype = genotyped.get_site()->get_genotype();
    GtypedIndices expected_gtype{2, 3};
    EXPECT_EQ(std::get<GtypedIndices>(gtype), expected_gtype);
}

TEST_F(TestLevelGenotyperModel_FourAlleles, GivenDifferentPloidies_CorrectNumberOfProducedGenotypes) {
    auto haploid_genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats);
    EXPECT_EQ(haploid_genotyped.get_likelihoods().size(), 4);

    auto diploid_genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats);
    // Expected number of genotypes: 4 diploid homozygous + (4 choose 2) diploid heterozygous
    EXPECT_EQ(diploid_genotyped.get_likelihoods().size(), 10);
}

/**
 * LevelGenotyper tests
 */
