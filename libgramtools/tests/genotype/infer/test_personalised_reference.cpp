#include "gtest/gtest.h"

#include "genotype/infer/output_specs/segment_tracker.hpp"
#include "genotype/infer/personalised_reference.hpp"
#include "prg/coverage_graph.hpp"

#include "../../test_resources/test_resources.hpp"
#include "mocks.hpp"

using namespace gram::genotype;

class Alleles_To_Paste : public ::testing::Test {
 protected:
  void SetUp() { site->set_alleles(all_alleles); }

  gt_site_ptr site = std::make_shared<MockGenotypedSite>();

  allele_vector all_alleles{
      Allele{"ATA", {0, 0, 0}, 0},
      Allele{"TTA", {0, 0, 0}, 1},
      Allele{"TTT", {0, 0, 0}, 2},
  };
};

TEST_F(Alleles_To_Paste, GivenInconsistentPloidy_Throws) {
  site->set_genotype(GtypedIndices{0, 1});
  EXPECT_THROW(get_all_alleles_to_paste(site, 3), InconsistentPloidyException);
}

TEST_F(Alleles_To_Paste, GivenGtype_CorrectAlleles) {
  site->set_genotype(GtypedIndices{0, 2});
  auto res = get_all_alleles_to_paste(site, 2);
  allele_vector expected{all_alleles.at(0), all_alleles.at(2)};
  EXPECT_EQ(res, expected);
}

TEST_F(Alleles_To_Paste, GivenNullGtype_CorrectAlleles) {
  site->set_genotype(GtypedIndices{-1});
  auto res = get_all_alleles_to_paste(site, 3);
  allele_vector expected{all_alleles.at(0), all_alleles.at(0),
                         all_alleles.at(0)};
  EXPECT_EQ(res, expected);
}

using str_vec = std::vector<std::string>;

class Personalised_Ref : public ::testing::Test {
 protected:
  void SetUp() {
    std::string linear_prg{"AT[CG[C,G]T,C]TT[AT,TT][C,G]"};
    PRG_String prg{prg_string_to_ints(linear_prg)};
    g = coverage_Graph{prg};
    graph_root = g.root;
    covG_ptrPair bubble;

    auto site1 = std::make_shared<MockGenotypedSite>();
    site1->set_alleles(allele_vector{
        Allele{"CGCT", {}, 0},
        Allele{"CGGT", {}, 0},
        Allele{"C", {}, 1},
    });
    bubble = get_bubble_nodes(g.bubble_map, 5);
    site1->set_site_end_node(bubble.second);

    // This, being nested, should get systematically skipped
    auto site2 = std::make_shared<MockGenotypedSite>();
    site2->set_alleles(allele_vector{
        Allele{"C", {}},
        Allele{"G", {}},
    });
    bubble = get_bubble_nodes(g.bubble_map, 7);
    site2->set_site_end_node(bubble.second);

    auto site3 = std::make_shared<MockGenotypedSite>();
    site3->set_alleles(allele_vector{
        Allele{"AT", {}},
        Allele{"TT", {}},
    });
    bubble = get_bubble_nodes(g.bubble_map, 9);
    site3->set_site_end_node(bubble.second);

    auto site4 = std::make_shared<MockGenotypedSite>();
    site4->set_alleles(allele_vector{
        Allele{"C", {}},
        Allele{"G", {}},
    });
    bubble = get_bubble_nodes(g.bubble_map, 11);
    site4->set_site_end_node(bubble.second);

    sites = gt_sites{site1, site2, site3, site4};
    set_trackers();
  }

  void null_all_sites() {
    // When all gts are null, ploidy is set to 1
    for (auto const& site : sites) site->set_genotype(GtypedIndices{-1});
  }

  void set_trackers() {
    std::stringstream ss{""};
    s1_tracker = SegmentTracker(ss);

    ss = std::stringstream{
        "chr1\t2\n"
        "chr2\t9\n"};
    s2_tracker_to_edge = SegmentTracker(ss);

    ss = std::stringstream{
        "chr1\t6\n"
        "chr2\t5\n"};
    s2_tracker_from_edge = SegmentTracker(ss);

    ss = std::stringstream{
        "chr1\t10\n"
        "chr2\t1\n"};
    s2_tracker_adjacentSites = SegmentTracker(ss);

    ss = std::stringstream{
        "chr1\t7\n"
        "chr2\t4\n"};
    s2_tracker_seq = SegmentTracker(ss);
  }

  SegmentTracker s1_tracker, s2_tracker_to_edge, s2_tracker_from_edge,
      s2_tracker_adjacentSites, s2_tracker_seq;
  coverage_Graph g;
  covG_ptr graph_root;
  gt_sites sites;
  str_vec res;
};

TEST_F(Personalised_Ref, GivenAllNullGts_CorrectInferredRef) {
  null_all_sites();
  auto result_map = get_personalised_ref(graph_root, sites, s1_tracker);
  auto result = *result_map.begin();
  std::string expected{"ATCGCTTTATC"};
  EXPECT_EQ(result.get_sequence(), expected);
}

TEST_F(Personalised_Ref, GivenHaploidGts_CorrectInferredRef) {
  sites.at(0)->set_genotype(GtypedIndices{2});
  sites.at(2)->set_genotype(GtypedIndices{1});
  sites.at(3)->set_genotype(GtypedIndices{1});
  auto result_map = get_personalised_ref(graph_root, sites, s1_tracker);
  auto result = *result_map.begin();
  std::string expected{"ATCTTTTG"};
  EXPECT_EQ(result.get_sequence(), expected);
}

TEST_F(Personalised_Ref, GivenHetDiploidGts_CorrectTwoInferredRefs) {
  sites.at(0)->set_genotype(GtypedIndices{1, 2});
  sites.at(2)->set_genotype(GtypedIndices{0, 1});
  sites.at(3)->set_genotype(GtypedIndices{0, 1});
  auto result_map = get_personalised_ref(graph_root, sites, s1_tracker);
  ASSERT_EQ(result_map.size(), 2);

  for (auto fasta : result_map) res.push_back(fasta.get_sequence());
  str_vec expected{{"ATCGGTTTATC"}, {"ATCTTTTG"}};
  EXPECT_EQ(res, expected);
}

TEST_F(Personalised_Ref, GivenHetSameGts_CorrectSingleInferredRef) {
  sites.at(0)->set_genotype(GtypedIndices{0, 0});
  sites.at(2)->set_genotype(GtypedIndices{1, 1});
  sites.at(3)->set_genotype(GtypedIndices{1, 1});
  auto result_vec = get_personalised_ref(graph_root, sites, s1_tracker);
  EXPECT_EQ(result_vec.size(), 2);

  unique_Fastas result_map{result_vec.begin(), result_vec.end()};

  auto first_ref = *result_map.begin();
  std::string expected_1{"ATCGCTTTTTG"};
  EXPECT_EQ(first_ref.get_sequence(), expected_1);
}

TEST_F(Personalised_Ref, GivenToEdgeS2Tracker_CorrectMultiSegRef) {
  null_all_sites();
  auto result_map = get_personalised_ref(graph_root, sites, s2_tracker_to_edge);
  ASSERT_EQ(result_map.size(), 2);

  for (auto fasta : result_map) res.push_back(fasta.get_sequence());
  str_vec expected{{"AT"}, {"CGCTTTATC"}};
  EXPECT_EQ(res, expected);
}

TEST_F(Personalised_Ref, GivenFromEdgeS2Tracker_CorrectMultiSegRef) {
  null_all_sites();
  auto result_map =
      get_personalised_ref(graph_root, sites, s2_tracker_from_edge);
  ASSERT_EQ(result_map.size(), 2);

  for (auto fasta : result_map) res.push_back(fasta.get_sequence());
  str_vec expected{{"ATCGCT"}, {"TTATC"}};
  EXPECT_EQ(res, expected);
}

TEST_F(Personalised_Ref, GivenAdjSitesS2Tracker_CorrectMultiSegRef) {
  null_all_sites();
  auto result_map =
      get_personalised_ref(graph_root, sites, s2_tracker_adjacentSites);
  ASSERT_EQ(result_map.size(), 2);

  for (auto fasta : result_map) res.push_back(fasta.get_sequence());
  str_vec expected{{"ATCGCTTTAT"}, {"C"}};
  EXPECT_EQ(res, expected);
}

TEST_F(Personalised_Ref, GivenSeqS2Tracker_CorrectMultiSegRef) {
  null_all_sites();
  auto result_map = get_personalised_ref(graph_root, sites, s2_tracker_seq);
  ASSERT_EQ(result_map.size(), 2);

  for (auto fasta : result_map) res.push_back(fasta.get_sequence());
  str_vec expected{{"ATCGCTT"}, {"TATC"}};
  EXPECT_EQ(res, expected);
}
