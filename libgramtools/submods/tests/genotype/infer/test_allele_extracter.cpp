#include "mocks.hpp"
#include "gtest/gtest.h"

#include "tests/common.hpp"
#include "genotype/infer/allele_extracter.hpp"
#include <memory>
#include <boost/shared_ptr.hpp>

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
    .WillRepeatedly(Return(GtypedIndices{0}));

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


TEST_F(AlleleCombineTest, TwoAllelesNullGenotype_oneCorrectCombinationAllele){

    EXPECT_CALL(site, get_genotype())
            .WillOnce(Return(false)); // Null genotype: should take the first allele in the vector

    EXPECT_CALL(site, get_alleles())
            .WillOnce(Return(
                    allele_vector{
                        Allele{ "TTT", {1, 1, 1} },
                        Allele{ "CCC", {0, 1, 1} }
                    }));

    allele_vector one_allele(existing_alleles.begin(), existing_alleles.begin() + 1);
    auto result = test_extracter.allele_combine(one_allele, 0);
    allele_vector expected {
            {
                    "ATTGTTT",
                    {0, 1, 2, 3, 1, 1, 1},
                    0
            }
    };

    EXPECT_EQ(result, expected);
};


TEST_F(AlleleCombineTest, TwoAllelesHeterozygousGenotype_FourCorrectCombinationAlleles){

    EXPECT_CALL(site, get_genotype())
            .WillRepeatedly(Return(GtypedIndices{0, 1}));

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


TEST(AllelePasteTest, TwoAllelesOneCoverageNode_CorrectlyAppendedSequenceandCoverage){

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

    // Note: need to explicitly pass in (dummy) site and allele IDs, else the Node thinks it is outside a variant site,
    // and does not need pb Coverage array.
    covG_ptr cov_Node = boost::make_shared<coverage_Node>("ATTCGC", 120, 1, 1);

    AlleleExtracter extracter;
    extracter.allele_paste(existing_alleles, cov_Node);

    allele_vector expected{
            {
                    "ATTGATTCGC",
                    {0, 1, 2, 3, 0, 0, 0, 0, 0, 0},
                    0
            },
            {
                    "ATCGATTCGC",
                    {0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
                    0
            }
    };

    EXPECT_EQ(existing_alleles, expected);
}

class AlleleExtracterTest : public ::testing::Test{
protected:
    std::shared_ptr<MockGenotypedSite> first_site_ptr = std::make_shared<MockGenotypedSite>();
    std::shared_ptr<MockGenotypedSite> second_site_ptr = std::make_shared<MockGenotypedSite>();
    gt_sites genotyped_sites = gt_sites{
            first_site_ptr,
            second_site_ptr
    };

    marker_vec v = prg_string_to_ints("AT[GCC[C,A,G]T,TTA]T");
    PRG_String prg_string{v};
    coverage_Graph cov_graph{prg_string};

    covG_ptrPair nested_bubble_nodes = get_bubble_nodes(cov_graph.bubble_map, 7);
    covG_ptrPair outer_bubble_nodes = get_bubble_nodes(cov_graph.bubble_map, 5);
};

TEST_F(AlleleExtracterTest, NestedBubble_CorrectAlleles){
    AlleleExtracter extracter{nested_bubble_nodes.first, nested_bubble_nodes.second, genotyped_sites};

    allele_vector expected{
            { "C", {0}, 0 },
            { "A", {0}, 1 },
            { "G", {0}, 2 }
    };
    EXPECT_EQ(extracter.get_alleles(), expected);
}

TEST_F(AlleleExtracterTest, OuterBubbleEncompassingHaploidNestedBubble_CorrectAlleles){
    EXPECT_CALL(*second_site_ptr, get_genotype())
    .WillOnce(Return(GtypedIndices{0}));

    EXPECT_CALL(*second_site_ptr, get_alleles())
    .WillOnce(Return(allele_vector{
            {"C", {0}, 0}
    }));

    EXPECT_CALL(*second_site_ptr, get_site_end_node())
    .WillOnce(Return(nested_bubble_nodes.second));

    AlleleExtracter extracter{outer_bubble_nodes.first, outer_bubble_nodes.second, genotyped_sites};

    allele_vector expected{
            { "GCCCT", {0, 0, 0, 0, 0}, 0 },
            { "TTA", {0, 0, 0}, 1 }
    };

    EXPECT_EQ(extracter.get_alleles(), expected);
}

TEST_F(AlleleExtracterTest, OuterBubbleEncompassingTriploidNestedBubble_CorrectAlleles){
    EXPECT_CALL(*second_site_ptr, get_genotype())
            .WillOnce(Return(GtypedIndices{0, 1, 2}));

    EXPECT_CALL(*second_site_ptr, get_alleles())
            .WillOnce(Return(allele_vector{
                    {"C", {0}, 0},
                    {"A", {0}, 1},
                    {"G", {0}, 2}
            }));

    EXPECT_CALL(*second_site_ptr, get_site_end_node())
            .WillOnce(Return(nested_bubble_nodes.second));

    AlleleExtracter extracter{outer_bubble_nodes.first, outer_bubble_nodes.second, genotyped_sites};

    allele_vector expected{
            { "GCCCT", {0, 0, 0, 0, 0}, 0 },
            { "GCCAT", {0, 0, 0, 0, 0}, 0 },
            { "GCCGT", {0, 0, 0, 0, 0}, 0 },
            { "TTA", {0, 0, 0}, 1 }
    };

    EXPECT_EQ(extracter.get_alleles(), expected);
}
