#include "gtest/gtest.h"
#include "prg/coverage_graph.hpp"
#include "prg/linearised_prg.hpp"
#include "simulate/induce_genotypes.hpp"
#include "simulate/simulate.hpp"
#include "test_resources/mocks.hpp"

using namespace gram::simulate;

class MakeRandomGenotypedSite : public ::testing::Test {
 protected:
  allele_vector alleles{
      Allele{"CTCGG", {}},
      Allele{"CG", {}},
      Allele{"CT", {}},
  };
  MockRandomGenerator rand;
};

using ::testing::Return;

TEST_F(MakeRandomGenotypedSite, GivenPickZerothAllele_CorrectSite) {
  EXPECT_CALL(rand, generate(0, 2)).WillOnce(Return(0));
  auto site = make_randomly_genotyped_site(&rand, alleles);
  allele_vector expected_als{alleles.begin(), alleles.begin() + 1};
  EXPECT_EQ(site->get_alleles(), expected_als);

  GtypedIndices expected_gts{0};
  EXPECT_EQ(site->get_genotype(), expected_gts);

  EXPECT_EQ(site->get_num_haplogroups(), 3);
}

TEST_F(MakeRandomGenotypedSite, GivenPickSecondAllele_CorrectSite) {
  EXPECT_CALL(rand, generate(0, 2)).WillOnce(Return(2));
  auto site = make_randomly_genotyped_site(&rand, alleles);
  allele_vector expected_als{alleles.at(0), alleles.at(2)};
  EXPECT_EQ(site->get_alleles(), expected_als);

  GtypedIndices expected_gts{1};  // Is rescaled
  EXPECT_EQ(site->get_genotype(), expected_gts);
}

TEST_F(MakeRandomGenotypedSite, GivenIgnoreREFAllele_CorrectSite) {
  EXPECT_CALL(rand, generate(1, 2)).WillOnce(Return(1));
  alleles.at(0).nesting_consistent = false;
  auto site = make_randomly_genotyped_site(&rand, alleles);
  allele_vector expected_als{alleles.at(0), alleles.at(1)};
  EXPECT_EQ(site->get_alleles(), expected_als);
}

class TestInduceGenotypes_ThreadSimpleSeq : public ::testing::Test {
 protected:
  coverage_Graph g;
  covG_ptr graph_end;
  TestInduceGenotypes_ThreadSimpleSeq() {
    auto encoded_prg = prg_string_to_ints("AA[A,C,G]TG[AC,[G,T]CA]CCC");
    PRG_String p{encoded_prg};
    g = coverage_Graph{p};
    auto cur_Node = g.root;
    while (cur_Node->get_num_edges() > 0)
      cur_Node = cur_Node->get_edges().at(0);
    graph_end = cur_Node;
  }
};

TEST_F(TestInduceGenotypes_ThreadSimpleSeq,
       GivenSequenceNotInGraph_ThrowsError) {
  std::string absent_sequence{"AACTGACTTT"};
  auto const result = thread_sequence(g.root, absent_sequence);
  EXPECT_EQ(result.size(), 0);
  EXPECT_THROW(get_single_endpoint(result, "", false), NoEndpoints);
}

TEST_F(TestInduceGenotypes_ThreadSimpleSeq,
       GivenSequenceInGraphButIncomplete_ThrowsError) {
  std::string sequence{"AACTGACC"};
  auto const result = thread_sequence(g.root, sequence);
  EXPECT_THROW(get_single_endpoint(result, "", false), NoEndpoints);
}

TEST_F(TestInduceGenotypes_ThreadSimpleSeq,
       GivenSequenceInGraphAndComplete_GetSingleEndpoint) {
  std::string goodseq1{"AACTGACCCC"};
  auto result = thread_sequence(g.root, goodseq1);
  EXPECT_EQ(result.size(), 1);
  auto endpoint = result.back();
  EXPECT_EQ(endpoint->get_prg_node(), graph_end);
  EXPECT_EQ(endpoint->get_offset(), 10);

  std::string goodseq2{"AAATGGCACCC"};
  result = thread_sequence(g.root, goodseq2);
  EXPECT_EQ(result.size(), 1);
  endpoint = result.back();
  EXPECT_EQ(endpoint->get_prg_node(), graph_end);
  EXPECT_EQ(endpoint->get_offset(), 11);
}

TEST(InduceGenotypes_ThreadAmbigSeq, FlexibleTreatmentOfAmbiguity) {
  // Below PRGs have sequence ambiguity
  auto ambiguous_prg = prg_string_to_ints("AA[A,AA]A[AA,A]");
  coverage_Graph g = coverage_Graph{PRG_String{ambiguous_prg}};

  auto endpoints = thread_sequence(g.root, "AAAAAA");
  EXPECT_TRUE(endpoints.size() > 1);
  EXPECT_THROW(get_single_endpoint(endpoints, "", true), TooManyEndpoints);

  // The last parameter to `get_single_endpoints` switches on/off tolerating
  // ambiguity
  ambiguous_prg = prg_string_to_ints("AT[CA,C[C,A]]GG");
  g = coverage_Graph{PRG_String{ambiguous_prg}};
  endpoints = thread_sequence(g.root, "ATCAGG");
  EXPECT_TRUE(endpoints.size() > 1);
  EXPECT_NO_THROW(get_single_endpoint(endpoints, "", false));
}

TEST(InduceGenotypes_NonConsumingInputSequence, LongestPathReturned) {
  // The threading process allows input sequences that consume the full graph
  // but not the full input sequence.
  // If there are several paths, return the most consuming one.
  std::vector<std::string> test_prg_seqs{"AA[A,AA]", "AA[AA,A]"};
  for (auto const& test_seq : test_prg_seqs) {
    auto linear_prg = prg_string_to_ints("AA[A,AA]");
    auto g = coverage_Graph{PRG_String{linear_prg}};
    auto endpoints = thread_sequence(g.root, "AAAAAAAA");
    EXPECT_EQ(endpoints.size(), 2);

    auto result = get_single_endpoint(endpoints, "", false);
    EXPECT_TRUE(result.first);  // `has_ambiguity` boolean
    EXPECT_EQ(result.second->get_offset(), 4);
  }
}

TEST(InduceGenotypes_ApplyGenotypes, GivenAmbiguousSequence_AMBIGFilterSet) {
  auto encoded_prg = prg_string_to_ints("AA[AA,A]A[A,AA]");
  auto g = coverage_Graph{PRG_String{encoded_prg}};
  auto sites = make_nulled_sites(g);
  auto endpoints = thread_sequence(g.root, "AAAAAA");
  auto checked_endpoint = get_single_endpoint(endpoints, "", false);
  apply_genotypes(checked_endpoint.second, checked_endpoint.first, sites);
  for (auto const& site : sites) EXPECT_TRUE(site->has_filter("AMBIG"));
}

TEST(InduceGenotypes_MakeNullSites, SitesAreNullGTAndHaveRefSeqOnly) {
  auto encoded_prg = prg_string_to_ints("AT[C,C[A,T]]GG");
  auto g = coverage_Graph{PRG_String{encoded_prg}};
  auto sites = make_nulled_sites(g);
  for (auto const& site : sites) {
    EXPECT_TRUE(site->is_null());
    EXPECT_EQ(site->get_alleles().size(), 1);
  }
  auto site1 = sites.at(0);
  EXPECT_EQ(site1->get_alleles().at(0).sequence, "C");

  auto site2 = sites.at(1);
  EXPECT_EQ(site2->get_alleles().at(0).sequence, "A");
}

class TestInduceGenotypes_InduceOneSeq : public ::testing::Test {
 protected:
  gt_sites sites;
  coverage_Graph g;
  void SetUp() {
    auto encoded_prg = prg_string_to_ints("AT[,C,GG]AA[TA,AA,G[GG,GGG]A,]CA");
    g = coverage_Graph{PRG_String{encoded_prg}};
    sites = make_nulled_sites(g);
  }
};

TEST_F(TestInduceGenotypes_InduceOneSeq,
       GivenRefThreadedSeq_CorrectGenotypedSites) {
  auto induced_sites = induce_genotypes_one_seq(sites, g, "ATAATACA", "");

  for (auto const& site :
       gt_sites{induced_sites.begin(), induced_sites.begin() + 2}) {
    EXPECT_FALSE(site->is_null());
    auto result = site->get_all_gtype_info();
    EXPECT_EQ(result.alleles.size(), 1);
    EXPECT_EQ(result.genotype, GtypedIndices{0});
    EXPECT_EQ(result.haplogroups, AlleleIds{0});
  }

  auto alleles = induced_sites.at(0)->get_alleles();
  EXPECT_EQ(alleles.at(0).sequence, "");

  alleles = induced_sites.at(1)->get_alleles();
  EXPECT_EQ(alleles.at(0).sequence, "TA");

  EXPECT_TRUE(induced_sites.at(2)->is_null());
}

TEST_F(TestInduceGenotypes_InduceOneSeq,
       GivenNonRefThreadedSeq_CorrectGenotypedSites) {
  auto induced_sites = induce_genotypes_one_seq(sites, g, "ATCAAGGGGACA", "");
  std::vector<std::string> expected_seqs;
  AlleleIds expected_ids;

  for (auto const& site : induced_sites) {
    EXPECT_FALSE(site->is_null());
    EXPECT_FALSE(site->has_filter("AMBIG"));
    auto result = site->get_all_gtype_info();
    EXPECT_EQ(result.alleles.size(), 2);
    EXPECT_EQ(result.genotype, GtypedIndices{1});
    EXPECT_EQ(result.haplogroups.size(), 1);

    expected_seqs.push_back(result.alleles.back().sequence);
    expected_ids.push_back(result.haplogroups.back());
  }

  EXPECT_EQ(expected_seqs, std::vector<std::string>({"C", "GGGGA", "GGG"}));
  EXPECT_EQ(expected_ids, AlleleIds({1, 2, 1}));
}
