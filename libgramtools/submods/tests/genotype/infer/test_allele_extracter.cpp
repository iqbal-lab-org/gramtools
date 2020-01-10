#include "mocks.hpp"
#include "gtest/gtest.h"

#include "genotype/infer/allele_extracter.hpp"
#include <memory>

using namespace ::testing;

class AlleleCombineTest : public ::testing::Test {
protected:

    std::shared_ptr<MockGenotypedSite> site_ptr = std::make_shared<MockGenotypedSite>();
    gt_sites sites{
            site_ptr
    };
    MockGenotypedSite& site = *site_ptr;
    AlleleExtracter test_extracter{sites};

    allele_vector existing_alleles{
            {
                    "ATTG",
                    {0, 1, 2, 3},
                    0
            },
            {
                    "ATCG",
                    {0, 0, 1, 1},
                    0
            }
    };
};

TEST_F(AlleleCombineTest, oneAlleleHaploidGenotype_oneCorrectCombinationAllele){

    EXPECT_CALL(site, get_genotype())
    .WillRepeatedly(Return(AlleleIds{0}));

    EXPECT_CALL(site, get_alleles())
    .WillRepeatedly(Return(
            allele_vector{Allele{
                    "CCC",
                    {1, 1, 1}
            }
    }));

    allele_vector one_allele(existing_alleles.begin(), existing_alleles.begin() + 1);
    auto result = test_extracter.allele_combine(one_allele, 0);
    allele_vector expected {
            {
                "ATTGCCC",
                {0, 1, 2, 3, 1, 1, 1},
                0
            }
    };

    EXPECT_EQ(result, expected);
};

TEST_F(AlleleCombineTest, TwoAllelesDiploidGenotype_FourCorrectCombinationAlleles){

    EXPECT_CALL(site, get_genotype())
            .WillRepeatedly(Return(AlleleIds{0, 1}));

    EXPECT_CALL(site, get_alleles())
            .WillRepeatedly(Return(
                    allele_vector{
                        Allele{
                            "CCC",
                            {1, 1, 1},
                            0
                        },
                        Allele{
                            "TTT",
                            {5, 5, 5},
                            1 // Note the pasted allele's haplogroup should get ignored
                        }
                    }));

    auto result = test_extracter.allele_combine(existing_alleles, 0);
    allele_vector expected {
            {
                    "ATTGCCC",
                    {0, 1, 2, 3, 1, 1, 1},
                    0
            },
            {
                    "ATTGTTT",
                    {0, 1, 2, 3, 5, 5, 5},
                    0
            },
            {
                    "ATCGCCC",
                    {0, 0, 1, 1, 1, 1, 1},
                    0
            },
            {
                    "ATCGTTT",
                    {0, 0, 1, 1, 5, 5, 5},
                    0
            },
    };

    EXPECT_EQ(result, expected);
};
