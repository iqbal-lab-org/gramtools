#include "gtest/gtest.h"
#include "mocks.hpp"
#include "genotype/infer/genotyped_site.hpp"

using ::testing::Return;

class GetUniqueGenotypedAlleles : public ::testing::Test{
protected:
    void SetUp(){
        EXPECT_CALL(site, get_alleles())
                .WillRepeatedly(Return( site_alleles ));
    }

    allele_vector site_alleles{
            Allele{"CCC", {1, 1, 1}},
            Allele{"GGG", {1, 1, 1}},
            Allele{"TTT", {1, 1, 1}},
    };
    MockGenotypedSite site;
};

TEST_F(GetUniqueGenotypedAlleles, GivenRepeatedGenotype_ProducedAllelesAreNotRepeated){
    EXPECT_CALL(site, get_genotype())
            .WillRepeatedly(Return(GtypedIndices{0, 0, 1}));

    auto extracted_alleles = site.extract_unique_genotyped_alleles();
    allele_vector expected{site_alleles.begin(), site_alleles.end() - 1};
    EXPECT_EQ(extracted_alleles, expected);
}

TEST_F(GetUniqueGenotypedAlleles, GivenNotUsingMock_ProducedAllelesSameAsWithMock){
    LevelGenotypedSite real_site;
    real_site.set_alleles(site_alleles);
    real_site.set_genotype(GtypedIndices{0, 0, 1}, 20);

    auto extracted_alleles = real_site.get_unique_genotyped_alleles();
    allele_vector expected{site_alleles.begin(), site_alleles.end() - 1};
    EXPECT_EQ(extracted_alleles, expected);
}

TEST_F(GetUniqueGenotypedAlleles, GivenUnorderedGenotype_ProducedAllelesAreOrdered) {
    EXPECT_CALL(site, get_genotype())
            .WillRepeatedly(Return(GtypedIndices{2, 0}));

    auto extracted_alleles = site.extract_unique_genotyped_alleles();
    allele_vector expected;
    expected.push_back(site_alleles.at(0));
    expected.push_back(site_alleles.at(2));
    EXPECT_EQ(extracted_alleles, expected);
}
