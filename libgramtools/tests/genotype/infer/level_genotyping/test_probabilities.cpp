#include "../mocks.hpp"
#include "genotype/infer/level_genotyping/runner.hpp"
#include "gtest/gtest.h"

using namespace gram::genotype::infer::probabilities;
using namespace ::testing;

TEST(ProbabilityMemoisation,
     GivenSameQueryParamsTwice_ProbabilityOnlyComputedOnce) {
  MockPmf pmf;

  EXPECT_CALL(pmf, compute_prob(params{1.5})).Times(1).WillOnce(Return(0.5));

  double prob1 = pmf(params{1.5});
  EXPECT_DOUBLE_EQ(prob1, 0.5);
  double prob2 = pmf(params{1.5});
  EXPECT_DOUBLE_EQ(prob2, prob1);
}

TEST(LikelihoodStats, DynamicChoiceOfProbDistribution) {
  // test based on failed dynamic casting producing nullptr, which evaluates to
  // false

  // here variance <= mean_cov_depth, so chooses Poisson
  auto lstats = LevelGenotyper::make_l_stats(10, 5, 0.01);
  EXPECT_FALSE(
      std::dynamic_pointer_cast<NegBinomLogPmf>(lstats.pmf_full_depth));
  EXPECT_TRUE(std::dynamic_pointer_cast<PoissonLogPmf>(lstats.pmf_full_depth));
  EXPECT_EQ(lstats.data_params, (DataParams{10., 0.01}));

  // here variance > mean cov depth, so chooses Negative Binomial
  lstats = LevelGenotyper::make_l_stats(10, 15, 0.01);
  EXPECT_FALSE(std::dynamic_pointer_cast<PoissonLogPmf>(lstats.pmf_full_depth));
  EXPECT_TRUE(std::dynamic_pointer_cast<NegBinomLogPmf>(lstats.pmf_full_depth));
}

TEST(LikelihoodStats, DynamicDataParams) {
  // tests lstats carries the parameters that are required for simulating data
  // according to chosen prob distributions

  // chooses Poisson: output mean_cov, pb_err_rate
  auto lstats = LevelGenotyper::make_l_stats(10, 5, 0.01);
  EXPECT_EQ(lstats.data_params, (DataParams{10., 0.01}));

  // chooses Negative Binomial: output mean_cov, num_successes, prob_success,
  // pb_err_rate
  lstats = LevelGenotyper::make_l_stats(10, 20, 0.01);
  DataParams expected(10., 0.01);
  expected.num_successes = 10.;
  expected.success_prob = 0.5;
  EXPECT_EQ(lstats.data_params, expected);

  // Test can get back original mean and variance from nbinom params
  auto num_successes = lstats.data_params.num_successes;
  auto prob_success = lstats.data_params.success_prob;
  EXPECT_EQ(int(num_successes * (1 - prob_success) / prob_success), 10);
  EXPECT_EQ(int(num_successes * (1 - prob_success) / pow(prob_success, 2)), 20);
}

TEST(LogPmfs, GivenConstructedObject_PmfAt0IsAlreadyMemoised) {
  pmf_ptr pmf;
  pmf = std::make_shared<PoissonLogPmf>(params{2});
  auto probs = pmf->get_probs();
  EXPECT_EQ(probs.size(), 1);
  EXPECT_EQ(probs.at(params{0}), -2);

  pmf = std::make_shared<NegBinomLogPmf>(params{2, 0.5});
  probs = pmf->get_probs();
  EXPECT_EQ(probs.size(), 1);
}

/*
 * truth probs computed using scipy 1.2.0 on Python 3.6.9.
 * Function: scipy.stats.<distrib>.pmf()
 */
TEST(LogPmfs, GivenTruthProbabilities_LogPmfValuesCorrect) {
  PoissonLogPmf dpois{params{2}};
  double known1{-1.3068528194400546};  // = ln(Poisson(lambda = 2, count = 2))
  auto res1 = dpois(params{2});
  EXPECT_FLOAT_EQ(res1, known1);

  dpois = PoissonLogPmf{params{2.5}};
  double known2{-1.3605657168116352};  // = ln(Poisson(lambda = 2, count = 2.5))
  auto res2 = dpois(params{2});
  EXPECT_DOUBLE_EQ(res2, known2);

  auto dnbinom = std::make_shared<NegBinomLogPmf>(params{2, 0.5});
  known1 = -1.6739764335716716;
  res1 = (*dnbinom)(params{2});
  EXPECT_DOUBLE_EQ(res1, known1);

  dnbinom = std::make_shared<NegBinomLogPmf>(params{2.5, 0.5});
  known2 = -2.3056313146033682;
  res2 = (*dnbinom)(params{4});
  EXPECT_DOUBLE_EQ(res2, known2);
}

TEST(MinCovMoreLikelyThanError,
     GivenMeanDepthAndErrorRate_CorrectMinCovThreshold) {
  LevelGenotyper g;
  std::vector<double> mean_depths{10, 10, 100};
  std::vector<double> pb_error_rates{0.0001, 0.001, 0.001};
  std::vector<CovCount> expected_min_cov_thresholds{1, 2, 10};

  std::size_t index{0};
  while (index < mean_depths.size()) {
    auto pmf = std::make_shared<PoissonLogPmf>(params{mean_depths.at(index)});
    auto min_cov_t =
        g.find_minimum_non_error_cov(pb_error_rates.at(index), pmf);
    EXPECT_EQ(min_cov_t, expected_min_cov_thresholds.at(index));
    ++index;
  }
}

TEST(MinCovMoreLikelyThanError, GivenZeroMeanDepth_ReturnsSomeValue) {
  LevelGenotyper g;
  // i.e. does not run into an infinite loop
  auto pmf = std::make_shared<PoissonLogPmf>(params{0});
  auto min_cov_t = g.find_minimum_non_error_cov(0.01, pmf);
  EXPECT_TRUE(min_cov_t >= 0);
}
