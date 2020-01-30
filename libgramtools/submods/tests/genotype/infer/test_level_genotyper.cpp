/**
 * Tests the internals of LevelGenotyperModel and of LevelGenotyper
 */
#include "gtest/gtest.h"
#include "mocks.hpp"
#include "genotype/infer/infer.hpp"

using namespace gram::genotype::infer;

TEST(HaploidCoverages, GivenSingletonCountsOnly_CorrectHaploidAndSingletonCovs){
   GroupedAlleleCounts  gp_covs{
                   {{0}, 5},
                   {{1}, 10},
                   {{3}, 1}
   };

   LevelGenotyperModel gtyper;
   gtyper.set_haploid_coverages(gp_covs, 4);
   PerAlleleCoverage expected_haploid_cov{5, 10, 0, 1};
   AlleleIdSet expected_singleton_cov{0, 1, 3};
   EXPECT_EQ(gtyper.get_haploid_covs(), expected_haploid_cov);
   EXPECT_EQ(gtyper.get_singleton_covs(), expected_singleton_cov);
}


TEST(HaploidCoverages, GivenMultiAllelicClasses_CorrectHaploidAndSingletonCovs){
    GroupedAlleleCounts  gp_covs{
            {{0}, 5},
            {{0, 1}, 4},
            {{1}, 10},
            {{2, 3}, 1}
    };

    LevelGenotyperModel gtyper;
    gtyper.set_haploid_coverages(gp_covs, 4);

    PerAlleleCoverage expected_haploid_cov{9, 14, 1, 1};
    AlleleIdSet expected_singleton_cov{0, 1};

    EXPECT_EQ(gtyper.get_haploid_covs(), expected_haploid_cov);
    EXPECT_EQ(gtyper.get_singleton_covs(), expected_singleton_cov);
}

TEST(DiploidCoverages, GivenMultiAllelicClasses_CorrectDiploidCovs){
    AlleleIds ids{0, 1}; // We want coverages of alleles 0 and 1

    GroupedAlleleCounts  gp_covs{
            {{0}, 7},
            {{0, 1}, 4},
            {{1}, 20},
            {{0, 3}, 3},
            {{2, 3}, 1}
    };

    // We have 10 units uniquely on 0, 20 uniquely on 1, and 4 shared between them.
    // These 4 should get dispatched in ratio 1:2 to alleles 0:1 (cf iqbal-lab-org/minos)

    LevelGenotyperModel gtyper;
    gtyper.set_haploid_coverages(gp_covs, 4);
    multiplicities haplogroup_multiplicities(4, false);
    auto diploid_covs = gtyper.compute_diploid_coverage(gp_covs, ids, haplogroup_multiplicities);
    EXPECT_FLOAT_EQ(diploid_covs.first, 10 + 4/3.);
    EXPECT_FLOAT_EQ(diploid_covs.second, 20 + 8/3.);
}

TEST(DiploidCoverages, GivenOnlyMultiAllelicClasses_CorrectDiploidCovs){
    AlleleIds ids{0, 1}; // We want coverages of alleles 0 and 1

    GroupedAlleleCounts  gp_covs{
            {{0, 1}, 3},
            {{2, 3}, 1}
    };

    // Edge case where singleton allele coverages are all 0
    // Then shared coverage should get dispatched equally (1:1 ratio)

    LevelGenotyperModel gtyper;
    gtyper.set_haploid_coverages(gp_covs, 4);
    multiplicities haplogroup_multiplicities(4, false);
    auto diploid_covs = gtyper.compute_diploid_coverage(gp_covs, ids, haplogroup_multiplicities);
    EXPECT_FLOAT_EQ(diploid_covs.first, 1.5);
    EXPECT_FLOAT_EQ(diploid_covs.second, 1.5);
}

class DiploidCoveragesOneDominatingClass : public ::testing::Test {
protected:
    void SetUp(){
        gtyper.set_haploid_coverages(gp_covs, 2);
    }

    GroupedAlleleCounts  gp_covs{
            {{0}, 8},
            {{0, 1}, 4},
    };

    LevelGenotyperModel gtyper;
};

TEST_F(DiploidCoveragesOneDominatingClass, GivenDifferentHaplogroups_CorrectDiploidCovs){
    // There is no unique coverage on haplogroup 1, thus all coverage goes to 0
    AlleleIds ids{0, 1};

    multiplicities haplogroup_multiplicities(2, false);
    auto diploid_covs = gtyper.compute_diploid_coverage(gp_covs, ids, haplogroup_multiplicities);
    EXPECT_FLOAT_EQ(diploid_covs.first, 12);
    EXPECT_FLOAT_EQ(diploid_covs.second, 0);
}

TEST_F(DiploidCoveragesOneDominatingClass, GivenSameHaplogroupTwice_CorrectDiploidCovs){
    // This can happen: when there is a nested site within, the extracted alleles have same haplogroup
    AlleleIds ids{0, 0};

    multiplicities haplogroup_multiplicities({true}); // The two alleles have the same haplogroup
    auto diploid_covs = gtyper.compute_diploid_coverage(gp_covs, ids, haplogroup_multiplicities);
    EXPECT_FLOAT_EQ(diploid_covs.first, 6);
    EXPECT_FLOAT_EQ(diploid_covs.second, 6);
}

TEST(CountCrediblePositions, GivenAlleleWithCredibleAndNonCrediblePositions_ReturnCrediblePositions){
   Allele test_allele{
      "ATCGCCG",
      {0, 0, 2, 3, 3, 5, 4},
      0
   };

   LevelGenotyperModel gtyper;
   auto num_credible = gtyper.count_credible_positions(3, test_allele);
   EXPECT_EQ(num_credible, 4);
}

TEST(CountTotalCov, GivenVariousCovStructures_CorrectTotalCoverages) {
    GroupedAlleleCounts gp_covs{};
    LevelGenotyperModel gtyper;
    EXPECT_EQ(gtyper.count_total_coverage(gp_covs), 0);

    GroupedAlleleCounts gp_covs2{
            {{0},    5},
            {{0, 1}, 4},
            {{1},    10},
            {{2, 3}, 1}
    };
    EXPECT_EQ(gtyper.count_total_coverage(gp_covs2), 20);
}

TEST(CountNumHaplogroups, GivenVariousAlleleVectors_CorrectNumHaplogroups){
    // Haplogroup should default to the same thing, consistently.
    allele_vector a1{
        Allele{"", {}},
        Allele{"", {}},
    };

    LevelGenotyperModel gtyper;
    multiplicities expected({true}); // Expect one entry, with more than one associated allele
    EXPECT_EQ(gtyper.count_num_haplogroups(a1), expected);

    allele_vector a2{
            Allele{"", {}, 0},
            Allele{"", {}, 1},
            Allele{"", {}, 1},
    };

    expected = multiplicities({false, true}); // Haplogroup 0 has 1 allele, haplogroup 1 has > 1 allele
    EXPECT_EQ(gtyper.count_num_haplogroups(a2), expected);
}

TEST(MakePermutations,GivenVariousParameters_CorrectPermutations){
    std::vector<GtypedIndices> expected;
    LevelGenotyperModel g;

    auto two_from_three = g.get_permutations(GtypedIndices{1,4,5}, 2);
    expected = {
            {1, 4},
            {1, 5},
            {4, 5}
    };
    EXPECT_EQ(two_from_three, expected);

    // Make sure result is internally sorted (at the genotype index level); needed for diploid coverage memoization
    auto from_unsorted = g.get_permutations(GtypedIndices{4,3,2}, 2);
    std::sort(from_unsorted.begin(), from_unsorted.end());
    expected = {
            {2, 3},
            {2, 4},
            {3, 4}
    };
    EXPECT_EQ(from_unsorted, expected);

    auto two_from_one = g.get_permutations(GtypedIndices{1}, 2); // Invalid call
    expected = {};
    EXPECT_EQ(two_from_one, expected);
}

TEST(RescaleGenotypes, GivenVariousGenotypes_CorrectRescaling){
    LevelGenotyperModel g;
    GtypedIndices no_zero_gt{1, 3};
    GtypedIndices no_zero_gt_rescaled{1, 2};

    EXPECT_EQ(g.rescale_genotypes(no_zero_gt), no_zero_gt_rescaled);

    GtypedIndices zero_and_repeated_gt{0, 4, 4};
    GtypedIndices  zero_and_repeated_gt_rescaled{0, 1, 1};
    EXPECT_EQ(g.rescale_genotypes( zero_and_repeated_gt ), zero_and_repeated_gt_rescaled);

    GtypedIndices shuffled_order{4, 2};
    GtypedIndices  shuffled_order_rescaled{1, 2};
    EXPECT_EQ(g.rescale_genotypes( shuffled_order ), shuffled_order_rescaled);
}


/**
 * Full run of the genotyping model
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
 * Test LevelGenotyper internals
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
