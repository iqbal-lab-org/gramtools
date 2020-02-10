/** @file
 *  - Tests genotyped site common interface routines
 */
#include "gtest/gtest.h"
#include "mocks.hpp"
#include "genotype/infer/level_genotyping/site.hpp"

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

TEST(NonGenotypedHaplogroups, GivenGenotypedSite_CorrectNonGenotypedHaplogroups){
    LevelGenotypedSite site;
    site.set_alleles(allele_vector{
        Allele{"ACGT", {1, 1, 1, 1}, 0},
        Allele{"TTTA", {1, 8, 1, 1}, 1},
        Allele{"TATA", {1, 8, 2, 1}, 1},
    });
    site.set_genotype(GtypedIndices{1, 2}, 5); // Het call of 2 alleles in same haplogroup.
    site.set_num_haplogroups(5);
    auto result = site.get_nonGenotyped_haplogroups();
    AlleleIds expected{0, 2, 3, 4};
    EXPECT_EQ(result, expected);
}

TEST(GetAllHaploGroups, GivenSiteWithGivenHaplotypeNum_CorrectReturnedHaplos){
    LevelGenotypedSite site;
    site.set_num_haplogroups(5);
    auto result = site.get_all_haplogroups();
    AlleleIds expected{0, 1, 2, 3, 4};
    EXPECT_EQ(result, expected);
}

TEST(GetGenotypedHaplogroups, GivenAllelesAndGT_CorrectHaplos){
    LevelGenotypedSite site;
    allele_vector alleles{
            Allele{"ACGT", {1, 1, 1, 1}, 0},
            Allele{"TTTA", {1, 8, 1, 1}, 1},
            Allele{"TATA", {1, 8, 2, 1}, 4},
    };
    GtypedIndices gt{0, 2};
    AlleleIds expected{0, 4};
    EXPECT_EQ(site.get_genotyped_haplogroups(alleles, gt), expected);
}