/**
 * Tests the internals of LevelGenotyperModel and of LevelGenotyper
 */
#include "genotype/infer/level_genotyping/model.hpp"
#include "genotype/infer/level_genotyping/runner.hpp"
#include "gtest/gtest.h"

using namespace gram::genotype::infer;

/**
 * Coverage computations
 */

TEST(HaploidCoverages,
     GivenSingletonCountsOnly_CorrectHaploidAndSingletonCovs) {
  GroupedAlleleCounts gp_covs{{{0}, 5}, {{1}, 10}, {{3}, 1}};

  LevelGenotyperModel gtyper;
  gtyper.set_haploid_coverages(gp_covs, 4);
  PerAlleleCoverage expected_haploid_cov{5, 10, 0, 1};
  EXPECT_EQ(gtyper.get_haploid_covs(), expected_haploid_cov);
  EXPECT_EQ(gtyper.get_singleton_covs(), expected_haploid_cov);
}

TEST(HaploidCoverages,
     GivenMultiAllelicClasses_CorrectHaploidAndSingletonCovs) {
  GroupedAlleleCounts gp_covs{{{0}, 5}, {{0, 1}, 4}, {{1}, 10}, {{2, 3}, 1}};

  LevelGenotyperModel gtyper;
  gtyper.set_haploid_coverages(gp_covs, 4);

  PerAlleleCoverage expected_haploid_cov{9, 14, 1, 1};
  PerAlleleCoverage expected_singleton_cov{5, 10, 0, 0};

  EXPECT_EQ(gtyper.get_haploid_covs(), expected_haploid_cov);
  EXPECT_EQ(gtyper.get_singleton_covs(), expected_singleton_cov);
}

TEST(DiploidCoverages, GivenMultiAllelicClasses_CorrectDiploidCovs) {
  AlleleIds ids{0, 1};  // We want coverages of alleles 0 and 1

  GroupedAlleleCounts gp_covs{
      {{0}, 7}, {{0, 1}, 4}, {{1}, 20}, {{0, 3}, 3}, {{2, 3}, 1}};

  // We have 10 units uniquely on 0, 20 uniquely on 1, and 4 shared between
  // them. These 4 should get dispatched in ratio 1:2 to alleles 0:1 (cf
  // iqbal-lab-org/minos)

  LevelGenotyperModel gtyper;
  gtyper.set_haploid_coverages(gp_covs, 4);
  multiplicities haplogroup_multiplicities(4, false);
  auto diploid_covs =
      gtyper.compute_diploid_coverage(gp_covs, ids, haplogroup_multiplicities);
  EXPECT_FLOAT_EQ(diploid_covs.first, 10 + 4 / 3.);
  EXPECT_FLOAT_EQ(diploid_covs.second, 20 + 8 / 3.);
}

TEST(DiploidCoverages, GivenOnlyMultiAllelicClasses_CorrectDiploidCovs) {
  AlleleIds ids{0, 1};  // We want coverages of alleles 0 and 1

  GroupedAlleleCounts gp_covs{{{0, 1}, 3}, {{2, 3}, 1}};

  // Edge case where singleton allele coverages are all 0
  // Then shared coverage should get dispatched equally (1:1 ratio)

  LevelGenotyperModel gtyper;
  gtyper.set_haploid_coverages(gp_covs, 4);
  multiplicities haplogroup_multiplicities(4, false);
  auto diploid_covs =
      gtyper.compute_diploid_coverage(gp_covs, ids, haplogroup_multiplicities);
  EXPECT_FLOAT_EQ(diploid_covs.first, 1.5);
  EXPECT_FLOAT_EQ(diploid_covs.second, 1.5);
}

TEST(LevelGenotyperModelDirectDeletion, GivenEmptyAllele_AssignsCoverage) {
  allele_vector alleles{
      Allele{"C", {8}, 0},
      Allele{"G", {8}, 0},
      Allele{"", {}, 1},
  };

  GroupedAlleleCounts gp_counts{
      {{0}, 8},
      {{1}, 8},
      {{0, 1}, 1},
  };
  auto expected = alleles;
  expected.at(2).pbCov = PerBaseCoverage{9};

  LevelGenotyperModel m;
  m.set_haploid_coverages(gp_counts, 2);
  m.assign_coverage_to_empty_alleles(alleles);
  EXPECT_EQ(alleles, expected);
}

class DiploidCoveragesOneDominatingClass : public ::testing::Test {
 protected:
  void SetUp() { gtyper.set_haploid_coverages(gp_covs, 2); }

  GroupedAlleleCounts gp_covs{
      {{0}, 8},
      {{0, 1}, 4},
  };

  LevelGenotyperModel gtyper;
};

TEST_F(DiploidCoveragesOneDominatingClass,
       GivenDifferentHaplogroups_CorrectDiploidCovs) {
  // There is no unique coverage on haplogroup 1, thus all coverage goes to 0
  AlleleIds ids{0, 1};

  multiplicities haplogroup_multiplicities(2, false);
  auto diploid_covs =
      gtyper.compute_diploid_coverage(gp_covs, ids, haplogroup_multiplicities);
  EXPECT_FLOAT_EQ(diploid_covs.first, 12);
  EXPECT_FLOAT_EQ(diploid_covs.second, 0);
}

TEST_F(DiploidCoveragesOneDominatingClass,
       GivenSameHaplogroupTwice_CorrectDiploidCovs) {
  // This can happen: when there is a nested site within, the extracted alleles
  // have same haplogroup
  AlleleIds ids{0, 0};

  multiplicities haplogroup_multiplicities(
      {true});  // The two alleles have the same haplogroup
  auto diploid_covs =
      gtyper.compute_diploid_coverage(gp_covs, ids, haplogroup_multiplicities);
  EXPECT_FLOAT_EQ(diploid_covs.first, 6);
  EXPECT_FLOAT_EQ(diploid_covs.second, 6);
}

TEST(CountCrediblePositions,
     GivenAlleleWithCredibleAndNonCrediblePositions_ReturnCredibleFraction) {
  Allele test_allele{"ATCGCCG", {0, 0, 2, 3, 3, 5, 4, 4}, 0};

  likelihood_related_stats l_stats;
  l_stats.credible_cov_t = 3;
  LevelGenotyperModel gtyper{l_stats, {}, {}};
  auto num_credible = gtyper.fraction_noncredible_positions(test_allele);
  EXPECT_DOUBLE_EQ(num_credible, 0.375);
}

TEST(CountTotalCov, GivenCovStructures_CorrectTotalCoverages) {
  GroupedAlleleCounts gp_covs{};
  LevelGenotyperModel gtyper;
  EXPECT_EQ(gtyper.count_total_coverage(gp_covs), 0);

  GroupedAlleleCounts gp_covs2{{{0}, 5}, {{0, 1}, 4}, {{1}, 10}, {{2, 3}, 1}};
  EXPECT_EQ(gtyper.count_total_coverage(gp_covs2), 20);
}

TEST(CountNumHaplogroups, GivenAlleleVectors_CorrectNumHaplogroups) {
  // Haplogroup should default to the same thing, consistently.
  allele_vector a1{
      Allele{"", {}},
      Allele{"", {}},
  };

  LevelGenotyperModel gtyper;
  multiplicities expected(
      {true});  // Expect one entry, with more than one associated allele
  EXPECT_EQ(gtyper.get_haplogroup_multiplicities(a1), expected);

  allele_vector a2{
      Allele{"", {}, 0},
      Allele{"", {}, 1},
      Allele{"", {}, 1},
  };

  expected = multiplicities(
      {false, true});  // Haplogroup 0 has 1 allele, haplogroup 1 has > 1 allele
  EXPECT_EQ(gtyper.get_haplogroup_multiplicities(a2), expected);
}

TEST(MakePermutations, GivenVariousParameters_CorrectPermutations) {
  std::vector<GtypedIndices> expected;
  LevelGenotyperModel g;

  auto two_from_three = g.get_permutations(GtypedIndices{1, 4, 5}, 2);
  expected = {{1, 4}, {1, 5}, {4, 5}};
  EXPECT_EQ(two_from_three, expected);

  // Make sure result is internally sorted (at the genotype index level); needed
  // for diploid coverage memoization
  auto from_unsorted = g.get_permutations(GtypedIndices{4, 3, 2}, 2);
  std::sort(from_unsorted.begin(), from_unsorted.end());
  expected = {{2, 3}, {2, 4}, {3, 4}};
  EXPECT_EQ(from_unsorted, expected);

  auto two_from_one = g.get_permutations(GtypedIndices{1}, 2);  // Invalid call
  expected = {};
  EXPECT_EQ(two_from_one, expected);
}

TEST(RescaleGenotypes, GivenVariousGenotypes_CorrectRescaling) {
  LevelGenotyperModel g;
  GtypedIndices no_zero_gt{1, 3};
  GtypedIndices no_zero_gt_rescaled{1, 2};

  EXPECT_EQ(g.rescale_genotypes(no_zero_gt), no_zero_gt_rescaled);

  GtypedIndices zero_and_repeated_gt{0, 4, 4};
  GtypedIndices zero_and_repeated_gt_rescaled{0, 1, 1};
  EXPECT_EQ(g.rescale_genotypes(zero_and_repeated_gt),
            zero_and_repeated_gt_rescaled);

  GtypedIndices shuffled_order{4, 2};
  GtypedIndices shuffled_order_rescaled{1, 2};
  EXPECT_EQ(g.rescale_genotypes(shuffled_order), shuffled_order_rescaled);
}

/**
 * Full run of the genotyping model
 */

TEST(TestLevelGenotyperModel_Failure, GivenOneAlleleOnly_Breaks) {
  // No likelihood ratio if only one allele. Note this should not present itself
  // if allele extraction works correctly, as any bubble has at least 2 alleles.
  allele_vector alleles{Allele{"ACGT", {1, 1, 1, 1}, 0}};
  GroupedAlleleCounts gp_counts;
  likelihood_related_stats l_stats;
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, false);
  EXPECT_DEATH(LevelGenotyperModel gtyper(data), "");
}

class TestLevelGenotyperModel_NullGTs : public ::testing::Test {
 protected:
  allele_vector alleles{
      Allele{"A", {0}, 0},
      Allele{"G", {0}, 1},
  };

  GroupedAlleleCounts gp_counts{};
  double mean_cov_depth{15}, mean_pb_error{0.01};
  likelihood_related_stats l_stats =
      LevelGenotyper::make_l_stats(mean_cov_depth, 0, mean_pb_error);
};

TEST_F(TestLevelGenotyperModel_NullGTs,
       GivenDuplicatedAllele_GTNullAndFilterSet) {
  alleles.emplace_back("A", PerBaseCoverage{1}, 1);
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, false);
  LevelGenotyperModel gtyper(data);
  EXPECT_TRUE(gtyper.get_site()->is_null());
  EXPECT_TRUE(gtyper.get_site()->has_filter("AMBIG"));
}

TEST_F(TestLevelGenotyperModel_NullGTs,
       Given0MeanCoverage_ReturnsNullGenotypedSite) {
  l_stats.data_params.mean_cov = 0;
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, false);
  auto genotyped = LevelGenotyperModel(data);

  EXPECT_TRUE(genotyped.get_site()->is_null());
}

TEST_F(TestLevelGenotyperModel_NullGTs,
       GivenNoCoverageOnAllAlleles_ReturnsNullGenotypedSite) {
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, false);
  auto genotyped = LevelGenotyperModel(data);

  EXPECT_TRUE(genotyped.get_site()->is_null());
}

TEST_F(TestLevelGenotyperModel_NullGTs,
       GivenSameCoverageOnAllAlleles_ReturnsNullGenotypedSite) {
  gp_counts = GroupedAlleleCounts{
      {{0}, 5},
      {{1}, 5},
  };
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, false);
  auto genotyped = LevelGenotyperModel(data);

  EXPECT_TRUE(genotyped.get_site()->is_null());
}

class TestLevelGenotyperModel_GTCalls : public ::testing::Test {
 protected:
  allele_vector alleles{
      Allele{"ATC", {0, 0, 1}, 0},
      Allele{"GGGCC", {10, 12, 12, 14, 14}, 1},
  };

  GroupedAlleleCounts gp_counts{
      {{0}, 1},
      {{1}, 13},
  };

  double mean_cov_depth{15}, mean_pb_error{0.01};
  likelihood_related_stats l_stats =
      LevelGenotyper::make_l_stats(mean_cov_depth, 0, mean_pb_error);
};

TEST_F(TestLevelGenotyperModel_GTCalls, CorrectDiploidCall) {
  ModelData data(alleles, gp_counts, Ploidy::Diploid, &l_stats, false);
  auto genotyped = LevelGenotyperModel(data);
  auto gtype = genotyped.get_site()->get_genotype();

  GtypedIndices expected_gtype{1, 1};
  EXPECT_EQ(gtype, expected_gtype);
}

TEST_F(TestLevelGenotyperModel_GTCalls, CorrectHaploidCall) {
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, false);
  auto genotyped = LevelGenotyperModel(data);
  auto gt_info = genotyped.get_site_gtype_info();
  EXPECT_EQ(gt_info.genotype, GtypedIndices{1});

  allele_vector expected_alleles{
      alleles.at(0),  // REF is not called, but still makes it in here
      alleles.at(1),
  };
  EXPECT_EQ(gt_info.alleles, expected_alleles);
}

TEST_F(TestLevelGenotyperModel_GTCalls,
       GivenLargeCovVariance_StillCorrectHaploidCall) {
  // Neg binomial gets used when variance cov depth exceed mean cov depth
  l_stats = LevelGenotyper::make_l_stats(mean_cov_depth, mean_cov_depth + 1,
                                         mean_pb_error);

  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, false);
  auto genotyped = LevelGenotyperModel(data);
  auto genotyped_nbinom = LevelGenotyperModel(data);
  auto gt_info = genotyped_nbinom.get_site_gtype_info();
  EXPECT_EQ(gt_info.genotype, GtypedIndices{1});
}

class TestLevelGenotyperModel_ExtraAlleles : public ::testing::Test {
 protected:
  allele_vector alleles{
      Allele{"A", {0}, 0},
      Allele{"G", {0}, 1},
  };
  // Allele at index 1 has highest log likelihood
  likelihood_map different_likelihoods_haploid{
      {-4, {0}},
      {-2, {1}},
  };
  likelihood_related_stats l_stats = LevelGenotyper::make_l_stats(40, 0, 0.01);
  multiplicities hap_muts{false, false};
};

TEST_F(TestLevelGenotyperModel_ExtraAlleles, GivenLargeCovs_NoExtraAlleles) {
  LevelGenotyperModel model(l_stats, {1, 39, 1}, different_likelihoods_haploid);
  model.CallGenotype(alleles, hap_muts, Ploidy::Haploid);
  EXPECT_FALSE(model.get_site()->extra_alleles());
}

TEST_F(TestLevelGenotyperModel_ExtraAlleles,
       Given0GtConf_AllBestAllelesAsExtraAlleles) {
  likelihood_map same_likelihoods{{-2, {0}}, {-2, {1}}};
  LevelGenotyperModel model(l_stats, {1, 39}, same_likelihoods);
  model.CallGenotype(alleles, hap_muts, Ploidy::Haploid);

  auto extra_alleles = model.get_site()->extra_alleles().value();
  EXPECT_EQ(extra_alleles, alleles);
  for (auto const& allele : extra_alleles)
    EXPECT_TRUE(allele.nesting_consistent);
}

TEST_F(TestLevelGenotyperModel_ExtraAlleles,
       GivenLowCovSituations_HaveExtraAlleles) {
  // Low total coverage on this site compared to mean cov of 40
  LevelGenotyperModel model(l_stats, {1, 5}, different_likelihoods_haploid);
  model.CallGenotype(alleles, hap_muts, Ploidy::Haploid);
  auto extra_alleles = model.get_site()->extra_alleles().value();
  EXPECT_EQ(extra_alleles, allele_vector{alleles.at(0)});
  EXPECT_FALSE(extra_alleles.at(0).nesting_consistent);

  // Low relative coverage difference on this site between alleles
  LevelGenotyperModel model2(l_stats, {20, 21}, different_likelihoods_haploid);
  model2.CallGenotype(alleles, hap_muts, Ploidy::Haploid);
  extra_alleles = model.get_site()->extra_alleles().value();
  EXPECT_EQ(extra_alleles, allele_vector{alleles.at(0)});
  EXPECT_FALSE(extra_alleles.at(0).nesting_consistent);
}

class TestLevelGenotyperModel_IgnoredREF : public ::testing::Test {
 public:
  void SetUp() {
    double mean_cov_depth{10}, mean_pb_error{0.01};
    l_stats = LevelGenotyper::make_l_stats(mean_cov_depth, 0, mean_pb_error);
  }

  allele_vector alleles{
      Allele{"A", {10}, 0, false},
      Allele{"C", {9}, 1},
      Allele{"G", {10}, 2},
  };

  GroupedAlleleCounts gp_counts{
      {{0}, 20},  // If ref considered, has more coverage
      {{1}, 9},
      {{2}, 10},
  };
  likelihood_related_stats l_stats;
};

TEST_F(TestLevelGenotyperModel_IgnoredREF, CorrectNumberOfLikelihoods) {
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, true);
  auto haploid_genotyped = LevelGenotyperModel(data);
  EXPECT_EQ(haploid_genotyped.get_likelihoods().size(), 2);

  data.ploidy = Ploidy::Diploid;
  auto diploid_genotyped = LevelGenotyperModel(data);
  // Two homs and one het
  EXPECT_EQ(diploid_genotyped.get_likelihoods().size(), 3);
}

TEST_F(TestLevelGenotyperModel_IgnoredREF,
       RescaledGenotypeIndices_AndREFInOutput) {
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, true);
  auto haploid_genotyped = LevelGenotyperModel(data);
  auto gtype_info = haploid_genotyped.get_site_gtype_info();
  allele_vector expected_alleles{alleles.at(0), alleles.at(2)};
  EXPECT_EQ(gtype_info.alleles, expected_alleles);
  GtypedIndices expected_gtypes{1};
  EXPECT_EQ(gtype_info.genotype, expected_gtypes);

  data.ploidy = Ploidy::Diploid;
  auto diploid_genotyped = LevelGenotyperModel(data);
  gtype_info = diploid_genotyped.get_site_gtype_info();
  // Expect heterozygous as homozygous calls have incompatible cov. counts
  EXPECT_EQ(gtype_info.alleles, alleles);
  expected_gtypes = GtypedIndices{1, 2};
  EXPECT_EQ(gtype_info.genotype, expected_gtypes);
}

TEST(TestLevelGenotyperModel, GivenHomozygousFavouredCov_GetHomozygousCall) {
  double mean_cov_depth{20}, mean_pb_error{0.01};

  allele_vector alleles{
      Allele{"AA", {0, 1}, 0},
      Allele{"TT", {20, 19}, 1},
  };

  GroupedAlleleCounts gp_counts{
      {{0}, 2},
      {{0, 1}, 1},
      {{1}, 20},
  };

  likelihood_related_stats l_stats =
      LevelGenotyper::make_l_stats(mean_cov_depth, 0, mean_pb_error);

  ModelData data(alleles, gp_counts, Ploidy::Diploid, &l_stats, false);
  auto genotyped = LevelGenotyperModel(data);
  auto gtype = genotyped.get_site()->get_genotype();
  GtypedIndices expected_gtype{1, 1};
  EXPECT_EQ(gtype, expected_gtype);
}

TEST(TestLevelGenotyperModel_NestedScenario,
     GivenAllelesInSameHaplogroup_CoverageGappedAllelePenalised) {
  allele_vector alleles{
      Allele{"AAAACAG", {0, 20, 20, 20, 20, 20, 0}, 0},
      Allele{"TAAACAT", {20, 20, 20, 20, 20, 20, 20}, 0},
  };
  GroupedAlleleCounts gp_counts{{{0}, 20}};

  // Deliberately enormous variance so it's plausible to have cov gaps
  // Nonetheless non-gapped allele is expected called due to gap penalty
  likelihood_related_stats l_stats =
      LevelGenotyper::make_l_stats(20, 200, 0.01);

  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats, false);
  auto genotyped = LevelGenotyperModel(data);
  EXPECT_EQ(genotyped.get_site()->get_genotype(), GtypedIndices{1});
}

TEST(TestLevelGenotyperModel_FourAlleles,
     GivenDifferentPloidies_CorrectNumberOfProducedGenotypes) {
  allele_vector alleles{
      Allele{"AATAA", {8, 8, 8, 8, 8}, 0}, Allele{"AAGAA", {7, 7, 7, 7, 7}, 0},
      Allele{"GGTGG",
             {15, 15, 15, 16, 16},
             1},  // 15 unique + 1 common with next allele
      Allele{"GGCGG",
             {14, 14, 14, 15, 15},
             1},  // 14 unique + 1 common with previous allele
  };

  GroupedAlleleCounts gp_counts{
      {{0}, 15},
      {{1}, 30},
  };
  double mean_cov_depth{30}, mean_pb_error{0.01};
  likelihood_related_stats l_stats =
      LevelGenotyper::make_l_stats(mean_cov_depth, 0, mean_pb_error);
  ModelData data(alleles, gp_counts, Ploidy::Haploid, &l_stats);
  auto haploid_genotyped = LevelGenotyperModel(data);
  EXPECT_EQ(haploid_genotyped.get_likelihoods().size(), 4);

  data.ploidy = Ploidy::Diploid;
  auto diploid_genotyped = LevelGenotyperModel(data);
  // Expected number of genotypes: 4 diploid homozygous + (4 choose 2) diploid
  // heterozygous
  EXPECT_EQ(diploid_genotyped.get_likelihoods().size(), 10);
}

class TestMaxLikelihoodCall : public ::testing::Test {
 public:
  likelihood_map likelihoods{
      {-1, {0}},
      {-2, {1}},
      {-3, {2}},
      {-4, {3}},
  };
  allele_vector alleles{
      Allele{"A", {}},
      Allele{"B", {}},
      Allele{"C", {}},
      Allele{"D", {}},
  };
};

TEST_F(TestMaxLikelihoodCall, LikelihoodsOrderedDescending) {
  auto it = likelihoods.begin();
  double max_likelihood = it->first;
  while (++it != likelihoods.end()) {
    EXPECT_TRUE(max_likelihood > it->first);
    max_likelihood = it->first;
  }
}

TEST_F(TestMaxLikelihoodCall, GivenOneLikelihood_Throws) {
  likelihood_map input{*likelihoods.begin()};
  EXPECT_THROW(LevelGenotyperModel::ChooseMaxLikelihood(input, allele_vector{}),
               IncorrectGenotyping);
}

TEST_F(TestMaxLikelihoodCall,
       GivenSeveralLikelihoods_ReturnsHighestLikelihood) {
  EXPECT_EQ(LevelGenotyperModel::ChooseMaxLikelihood(likelihoods, alleles),
            likelihoods.begin());
}

TEST_F(TestMaxLikelihoodCall, GivenInconsistentBestLikelihood_ItGetsSkipped) {
  alleles.at(0).nesting_consistent = false;
  auto expected = likelihoods.begin();
  ++expected;
  EXPECT_EQ(LevelGenotyperModel::ChooseMaxLikelihood(likelihoods, alleles),
            expected);
}

TEST_F(TestMaxLikelihoodCall,
       GivenInconsistentSecondBestLikelihood_NoSkipping) {
  alleles.at(1).nesting_consistent = false;
  EXPECT_EQ(LevelGenotyperModel::ChooseMaxLikelihood(likelihoods, alleles),
            likelihoods.begin());
}

TEST_F(TestMaxLikelihoodCall, GivenFewerThanTwoConsistentAlleles_Throws) {
  alleles.at(0).nesting_consistent = false;
  alleles.at(1).nesting_consistent = false;
  alleles.at(2).nesting_consistent = false;
  EXPECT_THROW(LevelGenotyperModel::ChooseMaxLikelihood(likelihoods, alleles),
               IncorrectGenotyping);
}

TEST_F(TestMaxLikelihoodCall, GivenNestingInconsistentBestAllele_NotCalled) {
  auto l_stats = LevelGenotyper::make_l_stats(20, 5, 0.01);
  multiplicities hap_muts{false};

  // Allele with most coverage is inconsistent
  // (this should not happen, but there are edge cases where it does: eg,
  // high-conf child SNP has per-base coverage closer to mean site coverage,
  // so it has higher likelihood)
  EXPECT_EQ(likelihoods.begin()->second, GtypedIndices{0});
  PerAlleleCoverage allele_covs{20, 15, 12, 8};
  alleles.at(0).nesting_consistent = false;

  LevelGenotyperModel model(l_stats, allele_covs, likelihoods);
  model.CallGenotype(alleles, hap_muts, Ploidy::Haploid);
  auto result = model.get_site_gtype_info();
  auto expected_alleles = allele_vector{alleles.at(0), alleles.at(1)};
  EXPECT_EQ(result.alleles, expected_alleles);
  EXPECT_EQ(result.genotype, GtypedIndices{1});
}
