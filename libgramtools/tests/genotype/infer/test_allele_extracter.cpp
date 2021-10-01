#include <memory>

#include "genotype/infer/allele_extracter.hpp"
#include "gtest/gtest.h"
#include "mocks.hpp"
#include "prg/coverage_graph.hpp"
#include "prg/types.hpp"
#include "submod_resources.hpp"

using namespace ::testing;
using namespace gram::submods;

TEST(ExtractRefAllele, GivenSiteNodesInGraph_CorrectRefAllele) {
  marker_vec v = prg_string_to_ints("AT[[C,A,G]T[G[,C]C,T],TTA]T");
  PRG_String prg_string{v};
  coverage_Graph cov_graph{prg_string};
  auto nodes = get_bubble_nodes(cov_graph.bubble_map, 5);
  auto ref_allele = extract_ref_allele(nodes.first, nodes.second);
  EXPECT_EQ(ref_allele.haplogroup, 0);
  EXPECT_EQ(ref_allele.sequence, "CTGC");
}

class AlleleCombineTest : public ::testing::Test {
 protected:
  AlleleCombineTest() {}

  std::shared_ptr<MockGenotypedSite> site_ptr =
      std::make_shared<MockGenotypedSite>();
  gt_sites sites{site_ptr};
  MockGenotypedSite& site = *site_ptr;
  AlleleExtracter test_extracter{sites};

  allele_vector existing_alleles{{"ATTG", {0, 1, 2, 3}, 0},
                                 {"ATCG", {0, 0, 1, 1}, 0}};
};

TEST_F(AlleleCombineTest,
       SiteHasOneCalledAllele_CorrectCombinationWithLeftHaplogroupKept) {
  // The called allele has hapg of 2, but we expect combined allele to keep its
  // haplogroup
  site.set_alleles(allele_vector{Allele{"CCC", {1, 1, 1}, 2}});
  site.set_genotype(GtypedIndices{0});

  allele_vector one_allele{existing_alleles.at(0)};
  auto result = test_extracter.allele_combine(one_allele, 0);
  allele_vector expected{{"ATTGCCC", {0, 1, 2, 3, 1, 1, 1}, 0}};
  EXPECT_EQ(result, expected);
}

TEST_F(AlleleCombineTest,
       SiteHasExtraAllele_ExtraAlleleIncludedAndNestingInconsistencyIncluded) {
  // Extraction includes extra alleles and nesting inconsistent gets copied to
  // combined allele
  site.set_alleles(allele_vector{
      Allele{"CCC", {1, 1, 1}},
      Allele{"GGG", {2, 2, 2}},
  });
  site.set_extra_alleles(allele_vector{Allele{"AAA", {2, 1, 0}, 2, false}});
  site.set_genotype(GtypedIndices{1});

  allele_vector one_allele{existing_alleles.at(0)};
  EXPECT_TRUE(one_allele.at(0).callable);

  auto result = test_extracter.allele_combine(one_allele, 0);
  allele_vector expected{
      {"ATTGGGG", {0, 1, 2, 3, 2, 2, 2}, 0},
      {"ATTGAAA", {0, 1, 2, 3, 2, 1, 0}, 0},
  };
  EXPECT_EQ(result, expected);
  EXPECT_TRUE(result.at(0).callable);
  EXPECT_FALSE(result.at(1).callable);
}

TEST_F(AlleleCombineTest, TwoAllelesNullGenotype_oneCorrectCombinationAllele) {
  site.set_genotype(GtypedIndices{-1});
  site.set_alleles(
      allele_vector{Allele{"TTT", {1, 1, 1}}, Allele{"CCC", {0, 1, 1}}});

  allele_vector one_allele(existing_alleles.begin(),
                           existing_alleles.begin() + 1);
  auto result = test_extracter.allele_combine(one_allele, 0);
  allele_vector expected{{"ATTGTTT", {0, 1, 2, 3, 1, 1, 1}, 0}};

  EXPECT_EQ(result, expected);
  EXPECT_TRUE(result.at(0).callable);
};

TEST_F(AlleleCombineTest,
       TwoAllelesHeterozygousGenotype_FourCorrectCombinationAlleles) {
  site.set_genotype(GtypedIndices{0, 1});

  site.set_alleles(allele_vector{
      Allele{"CCC", {1, 1, 1}, 0},
      Allele{
          "TTT",
          {5, 5, 5},
          1  // Note the pasted allele's haplogroup should get ignored
      }});

  auto result = test_extracter.allele_combine(existing_alleles, 0);
  allele_vector expected{
      {"ATTGCCC", {0, 1, 2, 3, 1, 1, 1}, 0},
      {"ATTGTTT", {0, 1, 2, 3, 5, 5, 5}, 0},
      {"ATCGCCC", {0, 0, 1, 1, 1, 1, 1}, 0},
      {"ATCGTTT", {0, 0, 1, 1, 5, 5, 5}, 0},
  };

  EXPECT_EQ(result, expected);
  for (auto const& allele : result) EXPECT_TRUE(allele.callable);
}

TEST(AllelePasteTest,
     TwoAllelesOneCoverageNode_CorrectlyAppendedSequenceandCoverage) {
  allele_vector existing_alleles{{"ATTG", {0, 1, 2, 3}, 0},
                                 {"ATCG", {0, 0, 1, 1}, 0}};

  // Note: need to explicitly pass in (dummy) site and allele IDs, else the Node
  // thinks it is outside a variant site, and does not need pb Coverage array.
  covG_ptr cov_Node = boost::make_shared<coverage_Node>("ATTCGC", 120, 1, 1);

  AlleleExtracter extracter;
  extracter.allele_paste(existing_alleles, cov_Node);

  allele_vector expected{{"ATTGATTCGC", {0, 1, 2, 3, 0, 0, 0, 0, 0, 0}, 0},
                         {"ATCGATTCGC", {0, 0, 1, 1, 0, 0, 0, 0, 0, 0}, 0}};

  EXPECT_EQ(existing_alleles, expected);
}

class AlleleExtracter_NestedPRG : public ::testing::Test {
 protected:
  void SetUp() {
    second_site_ptr->set_site_end_node(nested_bubble_nodes.second);
  }
  std::shared_ptr<MockGenotypedSite> first_site_ptr =
      std::make_shared<MockGenotypedSite>();
  std::shared_ptr<MockGenotypedSite> second_site_ptr =
      std::make_shared<MockGenotypedSite>();
  gt_sites genotyped_sites = gt_sites{first_site_ptr, second_site_ptr};

  marker_vec v = prg_string_to_ints("AT[GCC[C,A,G]T,TTA]T");
  PRG_String prg_string{v};
  coverage_Graph cov_graph{prg_string};

  covG_ptrPair nested_bubble_nodes = get_bubble_nodes(cov_graph.bubble_map, 7);
  covG_ptrPair outer_bubble_nodes = get_bubble_nodes(cov_graph.bubble_map, 5);
};

TEST_F(AlleleExtracter_NestedPRG, NestedBubble_CorrectAlleles) {
  AlleleExtracter extracter{nested_bubble_nodes.first,
                            nested_bubble_nodes.second, genotyped_sites};

  allele_vector expected{{"C", {0}, 0}, {"A", {0}, 1}, {"G", {0}, 2}};
  auto result = extracter.get_alleles();
  EXPECT_TRUE(result.at(0).callable);
  EXPECT_EQ(result, expected);
}

TEST_F(AlleleExtracter_NestedPRG,
       OuterBubbleEncompassingHaploidNestedBubble_CorrectAlleles) {
  second_site_ptr->set_genotype(GtypedIndices{0});
  second_site_ptr->set_alleles(allele_vector{{"C", {0}, 0}});

  AlleleExtracter extracter{outer_bubble_nodes.first, outer_bubble_nodes.second,
                            genotyped_sites};

  allele_vector expected{{"GCCCT", {0, 0, 0, 0, 0}, 0}, {"TTA", {0, 0, 0}, 1}};

  EXPECT_EQ(extracter.get_alleles(), expected);
}

TEST_F(AlleleExtracter_NestedPRG,
       OuterBubbleEncompassingTriploidNestedBubble_CorrectAlleles) {
  second_site_ptr->set_genotype(GtypedIndices{0, 1, 2});
  second_site_ptr->set_alleles(
      allele_vector{{"C", {0}, 0}, {"A", {0}, 1}, {"G", {0}, 2}});

  AlleleExtracter extracter{outer_bubble_nodes.first, outer_bubble_nodes.second,
                            genotyped_sites};

  allele_vector expected{{"GCCCT", {0, 0, 0, 0, 0}, 0},
                         {"GCCAT", {0, 0, 0, 0, 0}, 0},
                         {"GCCGT", {0, 0, 0, 0, 0}, 0},
                         {"TTA", {0, 0, 0}, 1}};

  auto result = extracter.get_alleles();
  EXPECT_TRUE(result.at(0).callable);
  EXPECT_EQ(result, expected);
}

TEST_F(AlleleExtracter_NestedPRG,
       OuterBubbleEncompassingHaploidNonREFNestedBubble_REFGetsProduced) {
  second_site_ptr->set_genotype(GtypedIndices{1});
  second_site_ptr->set_alleles(allele_vector{{"C", {0}, 0}, {"G", {0}, 2}});

  AlleleExtracter extracter{outer_bubble_nodes.first, outer_bubble_nodes.second,
                            genotyped_sites};

  // The REF (first allele in the site) needs to have gotten placed at index 0
  allele_vector expected{{"GCCCT", {0, 0, 0, 0, 0}, 0},
                         {"GCCGT", {0, 0, 0, 0, 0}, 0},
                         {"TTA", {0, 0, 0}, 1}};

  auto result = extracter.get_alleles();
  EXPECT_FALSE(result.at(0).callable);
  EXPECT_EQ(result, expected);
}

TEST_F(AlleleExtracter_NestedPRG,
       NestedBubbleHasNextBestAllele_NextBestAlleleGetsProduced) {
  second_site_ptr->set_genotype(GtypedIndices{1});
  second_site_ptr->set_alleles(allele_vector{{"C", {0}, 0}, {"G", {0}, 2}});
  second_site_ptr->set_extra_alleles(allele_vector{Allele{"A", {0}, 1}});

  AlleleExtracter extracter{outer_bubble_nodes.first, outer_bubble_nodes.second,
                            genotyped_sites};

  // The REF (first allele in the site) needs to have gotten placed at index 0
  allele_vector expected{{"GCCCT", {0, 0, 0, 0, 0}, 0},
                         {"GCCGT", {0, 0, 0, 0, 0}, 0},
                         {"GCCAT", {0, 0, 0, 0, 0}, 0},
                         {"TTA", {0, 0, 0}, 1}};

  EXPECT_EQ(extracter.get_alleles(), expected);
}

TEST(AlleleExtracter_DirectDeletionPRG,
     GivenOneBubble_DirectDeletionAlleleIsPresent) {
  marker_vec v = prg_string_to_ints("AT[GCC,TTA,]T");
  PRG_String prg_string{v};
  coverage_Graph cov_graph{prg_string};

  covG_ptrPair bubble_nodes = get_bubble_nodes(cov_graph.bubble_map, 5);
  gt_sites genotyped_sites;
  AlleleExtracter extracter{bubble_nodes.first, bubble_nodes.second,
                            genotyped_sites};

  allele_vector expected{
      Allele{"GCC", {0, 0, 0}, 0},
      Allele{"TTA", {0, 0, 0}, 1},
      Allele{"", {}, 2},
  };

  EXPECT_EQ(extracter.get_alleles(), expected);
}
