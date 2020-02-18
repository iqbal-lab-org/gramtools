#include "gtest/gtest.h"
#include "simulate/simulate.hpp"
#include "tests/common/mocks.hpp"

using namespace gram::simulate;

class MakeRandomGenotypedSite : public ::testing::Test{
protected:
    allele_vector alleles{
            Allele{"CTCGG", {}},
            Allele{"CG", {}},
            Allele{"CT", {}},
    };
    MockRandomGenerator rand;
};

using ::testing::Return;

TEST_F(MakeRandomGenotypedSite, GivenPickZerothAllele_CorrectSite){
   EXPECT_CALL(rand, generate(0, 2))
   .WillOnce(Return(0));
   auto site = make_randomly_genotyped_site(&rand, alleles);
   allele_vector expected_als{alleles.begin(), alleles.begin() + 1};
   EXPECT_EQ(site->get_alleles(), expected_als);

   GtypedIndices expected_gts{0};
   EXPECT_EQ(std::get<GtypedIndices>(site->get_genotype()), expected_gts);
}

TEST_F(MakeRandomGenotypedSite, GivenPickSecondAllele_CorrectSite){
    EXPECT_CALL(rand, generate(0, 2))
            .WillOnce(Return(2));
    auto site = make_randomly_genotyped_site(&rand, alleles);
    allele_vector expected_als{
        alleles.at(0),
        alleles.at(2)
    };
    EXPECT_EQ(site->get_alleles(), expected_als);

    GtypedIndices expected_gts{1}; // Is rescaled
    EXPECT_EQ(std::get<GtypedIndices>(site->get_genotype()), expected_gts);
}
