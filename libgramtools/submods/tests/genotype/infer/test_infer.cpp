/**
 * Tests full-round single-site genotyping:
 * - i. LevelGenotyper based one-site genotyping
 */
#include "gtest/gtest.h"
#include "genotype/infer/infer.hpp"
#include "mocks.hpp"

/**
 * LevelGenotyperModel tests (i.)
 */
TEST(TestLevelGenotyperModel_Failure, GivenOneAlleleOnly_Breaks){
    // No likelihood ratio if only one allele. Note this should not present itself if allele extraction
    // works correctly, as any bubble has at least 2 alleles.
    allele_vector alleles{
        Allele{"ACGT", {1, 1, 1, 1}, 0}
    };
    GroupedAlleleCounts gp_counts;
    likelihood_related_stats l_stats;
    EXPECT_DEATH(LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats, false), "");
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
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats, false);

    EXPECT_TRUE(genotyped.get_site()->is_null());
}

TEST_F(TestLevelGenotyperModel_NullGTs , GivenNoCoverageOnAllAlleles_ReturnsNullGenotypedSite) {
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats, false);

    EXPECT_TRUE(genotyped.get_site()->is_null());
}


class TestLevelGenotyperModel_TwoAllelesWithCoverage : public ::testing::Test {
protected:
    allele_vector alleles{
        Allele{ "ATCACC", {0, 0, 1, 1, 2, 2}, 0 },
        Allele{ "ATGACC", {0, 0, 0, 0, 1, 1}, 0 },
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
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats, false);

    auto genotyped_alleles = genotyped.get_site()->get_alleles();
    allele_vector expected_alleles{
            alleles.at(0), // REF is not called, but still makes it in here
            alleles.at(2),
    };
    EXPECT_EQ(genotyped_alleles, expected_alleles);

    auto gtype = genotyped.get_site()->get_genotype();
    EXPECT_TRUE(std::holds_alternative<GtypedIndices>(gtype));

    // The genotype needs to get rescaled: it is index 2 in original allele vector, but 1 in retained alleles
    GtypedIndices expected_gtype{1};
    EXPECT_EQ(std::get<GtypedIndices>(gtype), expected_gtype);
}


TEST_F(TestLevelGenotyperModel_TwoAllelesWithCoverage, GivenCoverage_ReturnsCorrectDiploidCall) {
    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats, false);

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

    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats, false);
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

    auto genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats, false);

    auto genotyped_alleles = genotyped.get_site()->get_alleles();
    allele_vector expected_alleles{
            alleles.at(0), // REF is not called, but still makes it in here
            alleles.at(2),
            alleles.at(3),
    };
    EXPECT_EQ(genotyped_alleles, expected_alleles);

    auto gtype = genotyped.get_site()->get_genotype();
    GtypedIndices expected_gtype{1, 2};
    EXPECT_EQ(std::get<GtypedIndices>(gtype), expected_gtype);
}

TEST_F(TestLevelGenotyperModel_FourAlleles, GivenDifferentPloidies_CorrectNumberOfProducedGenotypes) {
    auto haploid_genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats, false);
    EXPECT_EQ(haploid_genotyped.get_likelihoods().size(), 4);

    auto diploid_genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats, false);
    // Expected number of genotypes: 4 diploid homozygous + (4 choose 2) diploid heterozygous
    EXPECT_EQ(diploid_genotyped.get_likelihoods().size(), 10);
}

TEST(TestLevelGenotyperModel_IgnoredREF, GivenSeveralAllelesAndIgnoredREF_CorrectNumberOfLikelihoods){
    // Note this is only relevant to genotyping a site which contains a nested site
    // Here we imagine we called one allele in each of haplogroups 0 and 1, and that a REF got prepended
    allele_vector alleles{
            Allele{ "A", {0}, 0 },
            Allele{ "C", {8}, 0 },
            Allele{ "G", {8}, 1 },
    };

    GroupedAlleleCounts gp_counts{
            { {0}, 8 },
            { {1}, 8 },
    };
    double mean_cov_depth{8}, mean_pb_error{0.01};
    likelihood_related_stats l_stats = LevelGenotyper::make_l_stats(mean_cov_depth, mean_pb_error);

    // The last param to LevelGenotyperModel is whether to avoid using the REF
    auto haploid_genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Haploid, &l_stats, true);
    EXPECT_EQ(haploid_genotyped.get_likelihoods().size(), 2);

    auto diploid_genotyped = LevelGenotyperModel(&alleles, &gp_counts, Ploidy::Diploid, &l_stats, true);
    // Two homs and one het. Note this only works if you have singleton coverage on each haplogroup
    EXPECT_EQ(diploid_genotyped.get_likelihoods().size(), 3);
}


/**
 * LevelGenotyper tests (i.)
 */

TEST(LevelGenotyperInvalidation, GivenChildMapAndCandidateHaplos_CorrectHaplosWithSites){
    // site 7 lives on haplogroup 0 of site 5, and sites 9 and 11 live on its haplogroup 1.
    parental_map par_map{
            {7, VariantLocus{5, 1}},
            {9, VariantLocus{5, 2}},
            {11, VariantLocus{5, 2}},
    };
    child_map child_m = build_child_map(par_map);
    LevelGenotyper g{child_m, gt_sites{}};

    AlleleIds expected_haplogroups{0, 1}; // Expected in 0-based
    auto haplos_with_sites = g.get_haplogroups_with_sites(5, {0, 1, 2, 3});
    EXPECT_EQ(haplos_with_sites, expected_haplogroups);

    auto empty_query = g.get_haplogroups_with_sites(7, {0, 1, 2, 3});
    EXPECT_EQ(empty_query, AlleleIds{});
}

using ::testing::InSequence;
using ::testing::Return;
TEST(LevelGenotyperInvalidation, GivenNestingStructure_CorrectGenotypeNullifying){
    parental_map par_map{
            {7, VariantLocus{5, 1}},
            {9, VariantLocus{7, 2}},
    };
    child_map child_m = build_child_map(par_map);

    gt_sites sites(3);
    auto site1 = std::make_shared<MockGenotypedSite>();
    EXPECT_CALL(*site1, make_null())
            .Times(1);
    EXPECT_CALL(*site1, is_null())
            .WillOnce(Return(false));
    site1->set_num_haplogroups(5);
    sites.at(1) = site1;

    // SiteID 9 will get nulled by site 7.
    // Then when site 5 nulls site 7, I want site 9 to signal it is already nulled.
    auto site2 = std::make_shared<MockGenotypedSite>();
    {
        InSequence seq;
        EXPECT_CALL(*site2, is_null())
            .WillOnce(Return(false));
        EXPECT_CALL(*site2, make_null())
                .Times(1);
        EXPECT_CALL(*site2, is_null())
                .WillOnce(Return(true));
    }
    site2->set_num_haplogroups(5);
    sites.at(2) = site2;

    LevelGenotyper g{child_m, sites};

    // I want site 9 to be invalidated in this call
    g.invalidate_if_needed(7, AlleleIds{1});
    // And this call to invalidate site 7, without attempting to invalidate site 9 which was already so by above call.
    g.invalidate_if_needed(5, AlleleIds{0});
}
