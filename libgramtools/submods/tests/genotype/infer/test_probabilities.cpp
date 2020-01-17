#include "gtest/gtest.h"
#include "mocks.hpp"

using namespace gram::genotype::infer::probabilities;
using namespace ::testing;

TEST(ProbabilityMemoisation, GivenSameQueryParamsTwice_ProbabilityOnlyComputedOnce){
    MockPmf pmf;

    EXPECT_CALL(pmf, compute_prob(params{1.5}))
            .Times(1)
            .WillOnce(Return(0.5));

    float prob1 = pmf(params{1.5});
    EXPECT_FLOAT_EQ(prob1, 0.5);
    float prob2 = pmf(params{1.5});
    EXPECT_FLOAT_EQ(prob2, prob1);
}

TEST(PoissonLogPmf, GivenConstructedObject_Poisson0IsAlreadyMemoised){
    PoissonLogPmf pmf(params{2});
    auto probs = pmf.get_probs();
    EXPECT_EQ(probs.size(), 1);
    EXPECT_EQ(probs.at(params{0}), -2);
}

TEST(PoissonLogPmf, GivenKnownProbability_LogPoissonValueIsCorrect){
    /**
     * I computed the 'known' probs using Python 3.6.9, scipy 1.2.0. Function: scipy.stats.poisson.pmf()
     */
    PoissonLogPmf pmf(params{2});

    float known_log_poisson{-1.3068528194400546}; // = ln(Poisson(lambda = 2, count = 2))
    auto prob = pmf(params{2});
    EXPECT_FLOAT_EQ(prob, known_log_poisson);

    PoissonLogPmf pmf_float_mean(params{2.5});
    float known_log_poisson2{-1.3605657168116352}; // = ln(Poisson(lambda = 2, count = 2.5))
    auto prob2 = pmf_float_mean(params{2});
    EXPECT_FLOAT_EQ(prob2, known_log_poisson2);
}
