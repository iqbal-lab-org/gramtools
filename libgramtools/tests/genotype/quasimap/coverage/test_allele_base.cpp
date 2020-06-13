#include <cctype>

#include "gtest/gtest.h"

#include "../../../test_resources/test_resources.hpp"
#include "genotype/quasimap/coverage/allele_base.hpp"
#include "submod_resources.hpp"

using namespace gram::coverage::per_base;

TEST(AlleleBaseCoverageDump, GivenPopulatedAlleleBaseCoverage_CorrectJsonDump) {
  SitesAlleleBaseCoverage allele_base_coverage = {
      SitePbCoverage{
          PerBaseCoverage{1, 12},
          PerBaseCoverage{0, 3, 0},
      },
      SitePbCoverage{
          PerBaseCoverage{0},
          PerBaseCoverage{0, 19, 0},
      },
  };
  auto result = dump_allele_base_coverage(allele_base_coverage);
  std::string expected =
      "{\"allele_base_counts\":[[[1,12],[0,3,0]],[[0],[0,19,0]]]}";
  EXPECT_EQ(result, expected);
}

TEST(AlleleBaseCoverageDump,
     GivenSingleSiteAlleleBaseCoverage_CorrectJsonDump) {
  SitesAlleleBaseCoverage allele_base_coverage = {SitePbCoverage{
      PerBaseCoverage{1, 12},
      PerBaseCoverage{0, 3, 0},
  }};
  auto result = dump_allele_base_coverage(allele_base_coverage);
  std::string expected = "{\"allele_base_counts\":[[[1,12],[0,3,0]]]}";
  EXPECT_EQ(result, expected);
}

TEST(AlleleBaseCoverageDump, GivenEmptyAlleleBaseCoverage_CorrectJsonDump) {
  SitesAlleleBaseCoverage allele_base_coverage = {};
  auto result = dump_allele_base_coverage(allele_base_coverage);
  std::string expected = "{\"allele_base_counts\":[]}";
  EXPECT_EQ(result, expected);
}

TEST(AlleleBaseCoverageStructure, GivenNestedCovGraph_EmptyStructure) {
  auto prg_raw = prg_string_to_ints("[AC[TG,CC]T,T]A");
  auto prg_info = generate_prg_info(prg_raw);

  SitesAlleleBaseCoverage expected{};
  auto actual = coverage::generate::allele_base_non_nested(prg_info);
  EXPECT_EQ(actual, expected);
}

TEST(AlleleBaseCoverageStructure,
     GivenNonNestedCovGraphOneSite_CorrectStructure) {
  auto prg_raw = encode_prg("ac5gg6ga6ccc6c6aaa");
  auto prg_info = generate_prg_info(prg_raw);

  SitesAlleleBaseCoverage expected{
      SitePbCoverage{PerBaseCoverage{0, 0}, PerBaseCoverage{0, 0},
                     PerBaseCoverage{0, 0, 0}, PerBaseCoverage{0}}};
  auto actual = coverage::generate::allele_base_non_nested(prg_info);
  EXPECT_EQ(actual, expected);
}

TEST(AlleleBaseCoverageStructure,
     GivenNonNestedCovGraphTwoSitesAndOneEmptyAllele_CorrectStructure) {
  auto prg_raw = prg_string_to_ints("ac[a,c,tt]atg[gggg,,a]cc");
  auto prg_info = generate_prg_info(prg_raw);

  SitesAlleleBaseCoverage expected{SitePbCoverage{{0}, {0}, {0, 0}},
                                   SitePbCoverage{{0, 0, 0, 0}, {}, {0}}};
  auto actual = coverage::generate::allele_base_non_nested(prg_info);
  EXPECT_EQ(actual, expected);
}

TEST(DummyCovNode, BuildWithSizeSmallerThanEndCoord_ThrowsException) {
  EXPECT_THROW(DummyCovNode(0, 5, 3), InconsistentCovNodeCoordinates);
}

TEST(DummyCovNode, BuildWithStartGreaterThanEnd_ThrowsException) {
  EXPECT_THROW(DummyCovNode(2, 1, 3), InconsistentCovNodeCoordinates);
}

TEST(DummyCovNode, ExtendWithEndPosGreaterThanNodeSize_ThrowsException) {
  auto d = DummyCovNode(1, 5, 6);
  EXPECT_THROW(d.extend_coordinates({0, 6}), InconsistentCovNodeCoordinates);
}

TEST(DummyCovNode, ExtendNoStartAndNoEnd_CorrectUnchangesCoordinates) {
  auto d = DummyCovNode(1, 5, 6);
  d.extend_coordinates({2, 5});
  node_coordinates expected_coords{1, 5};
  EXPECT_EQ(expected_coords, d.get_coordinates());
}

TEST(DummyCovNode, ExtendStartAndEnd_CorrectExtendedCoordinates) {
  auto d = DummyCovNode(3, 3, 6);
  d.extend_coordinates({0, 5});
  node_coordinates expected_coords{0, 5};
  EXPECT_EQ(expected_coords, d.get_coordinates());
}

TEST(Traverser, StartOutOfSiteEndInSite_correctObjectState) {
  auto prg_raw = encode_prg("CT5gg6AAGa6cc");
  auto prg_info = generate_prg_info(prg_raw);

  std::size_t read_size = 5;
  VariantSitePath traversed_path{VariantLocus{5, FIRST_ALLELE + 1}};
  auto start_point = prg_info.coverage_graph.random_access[0];

  Traverser t{start_point, traversed_path, read_size};
  auto variant_node = t.next_Node().value();
  EXPECT_EQ(variant_node->get_site_ID(), 5);
  EXPECT_EQ(variant_node->get_allele_ID(), FIRST_ALLELE + 1);

  std::pair<uint32_t, uint32_t> expected_coordinates{0, 2};
  EXPECT_EQ(expected_coordinates, t.get_node_coordinates());
  EXPECT_EQ(false, t.next_Node().has_value());
}

TEST(Traverser, StartAndEndInSite_CorrectNodeInterval) {
  auto prg_raw = encode_prg("ct5g6aaAAAAAAaga6cc");
  auto prg_info = generate_prg_info(prg_raw);

  std::size_t read_size = 6;
  // Empty because the fact we are in VariantLocus{5, 2} is recorded in
  // traversing_path container
  VariantSitePath traversed_path{};
  auto start_point = prg_info.coverage_graph.random_access[7];

  Traverser t{start_point, traversed_path, read_size};
  auto variant_node = t.next_Node().value();

  std::pair<uint32_t, uint32_t> expected_coordinates{2, 7};
  EXPECT_EQ(expected_coordinates, t.get_node_coordinates());
}

TEST(Traverser, StartInSiteAndTraverseToAnotherSite_CorrectObjectState) {
  auto prg_raw = encode_prg("ct5g6aAA6CC7gc8ga8AAAAa8");
  auto prg_info = generate_prg_info(prg_raw);

  std::size_t read_size = 8;
  VariantSitePath traversed_path{VariantLocus{7, FIRST_ALLELE + 2}};
  auto start_point = prg_info.coverage_graph.random_access[6];

  Traverser t{start_point, traversed_path, read_size};
  auto cur_Node = t.next_Node();
  auto variant_node = cur_Node;
  while (cur_Node.has_value()) {
    variant_node = cur_Node;
    cur_Node = t.next_Node();
  }

  std::pair<uint32_t, uint32_t> expected_coordinates{0, 3};
  EXPECT_EQ(expected_coordinates, t.get_node_coordinates());
  EXPECT_EQ(0, t.get_remaining_bases());
}

// Helper function to get all the loci that were traversed. Modifies the
// traversal in place
VariantSitePath collect_traversal(Traverser& t) {
  VariantSitePath traversal;
  VariantLocus site_and_allele;
  auto cur_Node = t.next_Node();

  while (bool(cur_Node)) {
    site_and_allele = {cur_Node.value()->get_site_ID(),
                       cur_Node.value()->get_allele_ID()};
    traversal.push_back(site_and_allele);
    cur_Node = t.next_Node();
  }
  return traversal;
}

TEST(Traverser_Nested,
     StartOutOfSiteEndOutOfSite_CorrectChosenSitesAndEndState) {
  std::string raw_prg = "A[ctt,G[AAA,a]T]CCccc";
  marker_vec v = prg_string_to_ints(raw_prg);
  auto prg_info = generate_prg_info(v);

  std::size_t read_size = 8;
  VariantSitePath traversed_path{VariantLocus{7, FIRST_ALLELE},
                                 VariantLocus{5, FIRST_ALLELE + 1}};

  auto start_point = prg_info.coverage_graph.random_access[0];
  Traverser t{start_point, traversed_path, read_size};

  VariantSitePath expected_traversal{
      VariantLocus{5, FIRST_ALLELE + 1}, VariantLocus{7, FIRST_ALLELE},
      VariantLocus{
          5, FIRST_ALLELE + 1}  // After exiting site 7, we still have coverage
                                // to record on allele 2 of site 5 (base 'T')
  };

  VariantSitePath actual_traversal = collect_traversal(t);
  EXPECT_EQ(expected_traversal, actual_traversal);

  // Make sure we have consumed all bases of the read
  EXPECT_EQ(0, t.get_remaining_bases());
  // Make sure we are placed correctly in the last node
  std::pair<uint32_t, uint32_t> expected_last_node_coords{0, 1};
  EXPECT_EQ(expected_last_node_coords, t.get_node_coordinates());
}

TEST(Traverser_Nested,
     TraverseGraphWithLevel2Nesting_CorrectChosenSitesAndEndState) {
  std::string raw_prg = "A[CT[GC[c,A]A,gc]T[C,a]Tt,t]c";
  marker_vec v = prg_string_to_ints(raw_prg);
  auto prg_info = generate_prg_info(v);

  std::size_t read_size = 10;
  VariantSitePath traversed_path{
      VariantLocus{11, FIRST_ALLELE}, VariantLocus{9, FIRST_ALLELE + 1},
      VariantLocus{7, FIRST_ALLELE}, VariantLocus{5, FIRST_ALLELE}};
  auto start_point = prg_info.coverage_graph.random_access[0];
  Traverser t{start_point, traversed_path, read_size};

  VariantSitePath expected_traversal{
      VariantLocus{5, FIRST_ALLELE},     VariantLocus{7, FIRST_ALLELE},
      VariantLocus{9, FIRST_ALLELE + 1}, VariantLocus{7, FIRST_ALLELE},
      VariantLocus{5, FIRST_ALLELE},     VariantLocus{11, FIRST_ALLELE},
      VariantLocus{5, FIRST_ALLELE},
  };

  VariantSitePath actual_traversal = collect_traversal(t);
  EXPECT_EQ(expected_traversal, actual_traversal);

  EXPECT_EQ(0, t.get_remaining_bases());
  std::pair<uint32_t, uint32_t> expected_last_node_coords{0, 0};
  EXPECT_EQ(expected_last_node_coords, t.get_node_coordinates());
}

TEST(PbCovRecorder_NodeProcessing, ProcessNewCovNode_CorrectDummyCovNodeMade) {
  PbCovRecorder pb_rec;
  covG_ptr cov_node =
      boost::make_shared<coverage_Node>(coverage_Node{"ACTG", 102, 5, 2});
  realCov_to_dummyCov expected_mapping{{cov_node, DummyCovNode(1, 3, 4)}};

  pb_rec.process_Node(cov_node, 1, 3);
  EXPECT_EQ(expected_mapping, pb_rec.get_cov_mapping());
}

TEST(PbCovRecorder_NodeProcessing,
     ProcessExistingCovNode_CorrectlyUpdatedDummyCovNode) {
  covG_ptr cov_node =
      boost::make_shared<coverage_Node>(coverage_Node{"ACTGCC", 99, 5, 2});
  realCov_to_dummyCov existing_mapping{{cov_node, DummyCovNode{1, 3, 6}}};
  PbCovRecorder pb_rec(existing_mapping);
  pb_rec.process_Node(cov_node, 2, 5);

  realCov_to_dummyCov expected_mapping{{cov_node, DummyCovNode(1, 5, 6)}};

  EXPECT_EQ(expected_mapping, pb_rec.get_cov_mapping());
}

/**
 * Tests full coverage recording by inspecting `DummyCovNode`s and
 * `coverage_Node`s
 */
using dummy_cov_nodes = std::vector<DummyCovNode>;

dummy_cov_nodes collect_dummy_cov_nodes(coverage_Graph const& cov_graph,
                                        prg_positions positions,
                                        realCov_to_dummyCov cov_mapping) {
  dummy_cov_nodes result(positions.size());
  covG_ptr accessed_node;
  std::size_t index{0};

  for (auto& pos : positions) {
    accessed_node = cov_graph.random_access[pos].node;
    if (cov_mapping.find(accessed_node) == cov_mapping.end())
      result[index] = DummyCovNode{};
    else
      result[index] = cov_mapping.at(accessed_node);
    index++;
  }
  return result;
}

/*
PRG: GCT5C6G6T6AG7T8CC8CT
i	BWT	SA	text_suffix
0	T	20
1	6	10	A G 7 T 8 C C 8 C T
2	8	15	C C 8 C T
3	8	18	C T
4	G	1	C T 5 C 6 G 6 T 6 A G 7 T 8 C C 8 C T
5	5	4	C 6 G 6 T 6 A G 7 T 8 C C 8 C T
6	C	16	C 8 C T
7	0	0	G C T 5 C 6 G 6 T 6 A G 7 T 8 C C 8 C T
8	6	6	G 6 T 6 A G 7 T 8 C C 8 C T
9	A	11	G 7 T 8 C C 8 C T
10	C	19	T
11	C	2	T 5 C 6 G 6 T 6 A G 7 T 8 C C 8 C T
12	6	8	T 6 A G 7 T 8 C C 8 C T
13	7	13	T 8 C C 8 C T
14	T	3	5 C 6 G 6 T 6 A G 7 T 8 C C 8 C T
15	T	9	6 A G 7 T 8 C C 8 C T
16	C	5	6 G 6 T 6 A G 7 T 8 C C 8 C T
17	G	7	6 T 6 A G 7 T 8 C C 8 C T
18	G	12	7 T 8 C C 8 C T
19	T	14	8 C C 8 C T
20	C	17	8 C T
*/
class PbCovRecorder_TwoSitesNoNesting : public ::testing::Test {
 protected:
  void SetUp() {
    std::string raw_prg = "GCT5C6G6T6AG7T8CC8CT";
    marker_vec v = encode_prg(raw_prg);
    prg_info = generate_prg_info(v);
  }

  PRG_Info prg_info;
  prg_positions all_sequence_node_positions{0, 4, 6, 8, 10, 13, 15, 18};

  // Read: CTGAGC from pos 1
  std::size_t read1_size = 6;
  SearchState read_1{SA_Interval{4, 4}, VariantSitePath{
                                            VariantLocus{7, FIRST_ALLELE + 1},
                                            VariantLocus{5, FIRST_ALLELE + 1},
                                        }};

  // Read: TAGCCC from pos 8
  std::size_t read2_size = 7;
  SearchState read_2{SA_Interval{12, 12}, VariantSitePath{
                                              VariantLocus{7, FIRST_ALLELE + 1},
                                          }};
};

TEST_F(PbCovRecorder_TwoSitesNoNesting,
       ReadCoversTwoSites_CorrectCoverageNodes) {
  // PRG: "gCT5c6G6t6AG7t8Cc8ct" ; Read: "CTGAGC"
  PbCovRecorder{prg_info, SearchStates{read_1}, read1_size};
  auto actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

  SitePbCoverage expected_coverage{PerBaseCoverage{},     PerBaseCoverage{0},
                                   PerBaseCoverage{1},    PerBaseCoverage{0},
                                   PerBaseCoverage{},     PerBaseCoverage{0},
                                   PerBaseCoverage{1, 0}, PerBaseCoverage{}};

  EXPECT_EQ(expected_coverage, actual_coverage);
}

TEST_F(PbCovRecorder_TwoSitesNoNesting,
       ReadCoversTwoSites2_CorrectCoverageNodes) {
  // PRG: "GCT5C6G6T6AG7T8CC8CT" ; Read: "TAGCCC"

  PbCovRecorder{prg_info, SearchStates{read_2}, read2_size};
  auto actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

  SitePbCoverage expected_coverage{PerBaseCoverage{},     PerBaseCoverage{0},
                                   PerBaseCoverage{0},    PerBaseCoverage{1},
                                   PerBaseCoverage{},     PerBaseCoverage{0},
                                   PerBaseCoverage{1, 1}, PerBaseCoverage{}};

  EXPECT_EQ(expected_coverage, actual_coverage);
}

/*
PRG: AAT[ATAT,AA,]AGG
i	BWT	SA	text_suffix
0	G	16	0
1	0	0	A A T 5 A T A T 6 A A 6 6 A G G 0
2	6	9	A A 6 6 A G G 0
3	6	13	A G G 0
4	5	4	A T A T 6 A A 6 6 A G G 0
5	A	1	A T 5 A T A T 6 A A 6 6 A G G 0
6	T	6	A T 6 A A 6 6 A G G 0
7	A	10	A 6 6 A G G 0
8	G	15	G 0
9	A	14	G G 0
10	A	5	T A T 6 A A 6 6 A G G 0
11	A	2	T 5 A T A T 6 A A 6 6 A G G 0
12	A	7	T 6 A A 6 6 A G G 0
13	T	3	5 A T A T 6 A A 6 6 A G G 0
14	T	8	6 A A 6 6 A G G 0
15	6	12	6 A G G 0
16	A	11	6 6 A G G 0
 */
class PbCovRecorder_WithRepeatsAndEmptyAllele : public ::testing::Test {
 protected:
  void SetUp() {
    std::string raw_prg = "AAT[ATAT,AA,]AGG";
    marker_vec v = prg_string_to_ints(raw_prg);
    prg_info = generate_prg_info(v);
  }
  PRG_Info prg_info;
  prg_positions all_sequence_node_positions{0, 4, 9, 12};

  // Read: ATAT, occurs twice: from pos 1 and from pos 4
  std::size_t read1_size = 4;
  SearchStates read_1{
      SearchState{SA_Interval{4, 4}, VariantSitePath{}},
      SearchState{SA_Interval{5, 5},
                  VariantSitePath{VariantLocus{5, FIRST_ALLELE}}}};

  // Read: ATAAA, occurs from pos 1
  std::size_t read2_size = 5;
  SearchState read_2{SearchState{
      SA_Interval{5, 5}, VariantSitePath{VariantLocus{5, FIRST_ALLELE + 1}}}};

  // Read: AATAGG, occurs from pos 0, goes through deletion
  std::size_t read3_size = 5;
  SearchState read_3{SearchState{
      SA_Interval{1, 1}, VariantSitePath{VariantLocus{5, FIRST_ALLELE + 2}}}};
};

TEST_F(PbCovRecorder_WithRepeatsAndEmptyAllele,
       RepeatedMultiMappedRead_CoverageOnlyAddedOnce) {
  // PRG: "AAT[ATAT,AA,]AGG" ; Read: ATAT
  PbCovRecorder{prg_info, read_1, read1_size};
  auto actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

  SitePbCoverage expected_coverage{PerBaseCoverage{},
                                   PerBaseCoverage{1, 1, 1, 1},
                                   PerBaseCoverage{0, 0}, PerBaseCoverage{}};

  EXPECT_EQ(expected_coverage, actual_coverage);
}

TEST_F(PbCovRecorder_WithRepeatsAndEmptyAllele,
       MapAReadMultipleSeparateTimes_CoverageCorrectlyMultiplyAdded) {
  // PRG: "AAT[ATAT,AA]AGG" ; Read: ATAAA
  uint16_t i;
  for (i = 0; i <= 2; i++)
    PbCovRecorder{prg_info, SearchStates{read_2}, read2_size};
  auto actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

  SitePbCoverage expected_coverage{PerBaseCoverage{},
                                   PerBaseCoverage{0, 0, 0, 0},
                                   PerBaseCoverage{3, 3}, PerBaseCoverage{}};

  EXPECT_EQ(expected_coverage, actual_coverage);

  // Collect coverage on the deletion read: ATAGG
  // No pb coverage recorded for it as it is not represented as a node
  for (i = 0; i <= 4; i++)
    PbCovRecorder{prg_info, SearchStates{read_3}, read3_size};
  actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);
  EXPECT_EQ(expected_coverage, actual_coverage);
}

/*
PRG: AT[GC[GCC,CCGC],T]TTTT
i	BWT	SA	text_suffix
0	T	22	0
1	0	0	A T 5 G C 7 G C C 8 C C G C 8 6 T 6 T T T T 0
2	8	10	C C G C 8 6 T 6 T T T T 0
3	G	7	C C 8 C C G C 8 6 T 6 T T T T 0
4	C	11	C G C 8 6 T 6 T T T T 0
5	G	4	C 7 G C C 8 C C G C 8 6 T 6 T T T T 0
6	C	8	C 8 C C G C 8 6 T 6 T T T T 0
7	G	13	C 8 6 T 6 T T T T 0
8	7	6	G C C 8 C C G C 8 6 T 6 T T T T 0
9	5	3	G C 7 G C C 8 C C G C 8 6 T 6 T T T T 0
10	C	12	G C 8 6 T 6 T T T T 0
11	T	21	T 0
12	T	20	T T 0
13	T	19	T T T 0
14	6	18	T T T T 0
15	A	1	T 5 G C 7 G C C 8 C C G C 8 6 T 6 T T T T 0
16	6	16	T 6 T T T T 0
17	T	2	5 G C 7 G C C 8 C C G C 8 6 T 6 T T T T 0
18	T	17	6 T T T T 0
19	8	15	6 T 6 T T T T 0
20	C	5	7 G C C 8 C C G C 8 6 T 6 T T T T 0
21	C	9	8 C C G C 8 6 T 6 T T T T 0
22	C	14	8 6 T 6 T T T T 0
*/
class PbCovRecorder_nestedDeletion : public ::testing::Test {
 protected:
  void SetUp() {
    std::string raw_prg = "AT[GC[GCC,CCGC],T]TTTT";
    marker_vec v = prg_string_to_ints(raw_prg);
    prg_info = generate_prg_info(v);
  }
  PRG_Info prg_info;
  prg_positions all_sequence_node_positions{0, 3, 6, 10, 16, 18};

  // Make some read SearchStates
  // Read: CGCCTT
  SearchState simple_read_1{SA_Interval{5, 5},
                            VariantSitePath{VariantLocus{7, FIRST_ALLELE}}};

  // Read: ATTTT
  SearchState simple_read_2{SA_Interval{1, 1},
                            VariantSitePath{VariantLocus{5, FIRST_ALLELE + 1}}};

  // Read: GCC. Two distinct occurrences compatible with same sites
  SearchStates multi_mapped_reads_1{
      SearchState{SA_Interval{9, 9},
                  VariantSitePath{VariantLocus{7, FIRST_ALLELE + 1}}},
      SearchState{SA_Interval{8, 8}, VariantSitePath{}}};

  // Read: CTT. In a single Search State
  SearchStates multi_mapped_reads_2{
      SearchState{SA_Interval{6, 7}, VariantSitePath{}}};
};

TEST_F(PbCovRecorder_nestedDeletion, simpleRead1Mapped_correctDummyCovNodes) {
  // PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "CGCCTT"
  std::size_t read_size{6};
  PbCovRecorder recorder(prg_info, read_size);
  recorder.process_SearchState(simple_read_1);
  auto actual_dummies = collect_dummy_cov_nodes(prg_info.coverage_graph,
                                                all_sequence_node_positions,
                                                recorder.get_cov_mapping());

  dummy_cov_nodes expected_dummies{DummyCovNode{},        DummyCovNode{1, 1, 2},
                                   DummyCovNode{0, 2, 3}, DummyCovNode{},
                                   DummyCovNode{},        DummyCovNode{}};
  EXPECT_EQ(expected_dummies, actual_dummies);
}

TEST_F(PbCovRecorder_nestedDeletion,
       simpleRead1Mapped_correctRecordedPbCoverage) {
  // PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "CGCCTT"
  SearchStates mapping{simple_read_1};
  std::size_t read_size{6};
  PbCovRecorder recorder(prg_info, mapping, read_size);
  auto actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

  SitePbCoverage expected_coverage{
      PerBaseCoverage{},        PerBaseCoverage{0, 1},
      PerBaseCoverage{1, 1, 1}, PerBaseCoverage{0, 0, 0, 0},
      PerBaseCoverage{0},       PerBaseCoverage{}};
  EXPECT_EQ(expected_coverage, actual_coverage);
}

TEST_F(PbCovRecorder_nestedDeletion, simpleRead2Mapped_correctDummyCovNodes) {
  // PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "ATTTT"
  SearchStates mapping = SearchStates{simple_read_2};
  std::size_t read_size{5};
  PbCovRecorder recorder(prg_info, read_size);
  recorder.process_SearchState(simple_read_2);
  auto actual_dummies = collect_dummy_cov_nodes(prg_info.coverage_graph,
                                                all_sequence_node_positions,
                                                recorder.get_cov_mapping());

  dummy_cov_nodes expected_dummies{DummyCovNode{},        DummyCovNode{},
                                   DummyCovNode{},        DummyCovNode{},
                                   DummyCovNode{0, 0, 1}, DummyCovNode{}};
  EXPECT_EQ(expected_dummies, actual_dummies);
}

TEST_F(PbCovRecorder_nestedDeletion,
       simpleRead2Mapped_correctRecordedPbCoverage) {
  // PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "ATTTT"
  SearchStates mapping{simple_read_2};
  std::size_t read_size{5};
  PbCovRecorder recorder(prg_info, mapping, read_size);
  auto actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

  SitePbCoverage expected_coverage{
      PerBaseCoverage{},        PerBaseCoverage{0, 0},
      PerBaseCoverage{0, 0, 0}, PerBaseCoverage{0, 0, 0, 0},
      PerBaseCoverage{1},       PerBaseCoverage{}};
  EXPECT_EQ(expected_coverage, actual_coverage);
}

TEST_F(PbCovRecorder_nestedDeletion,
       multiMappedReadDistinctSearchStates_correctRecordedPbCoverage) {
  // PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "GCC"
  std::size_t read_size{3};
  PbCovRecorder{prg_info, multi_mapped_reads_1, read_size};
  auto actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

  SitePbCoverage expected_coverage{
      PerBaseCoverage{},        PerBaseCoverage{1, 1},
      PerBaseCoverage{1, 1, 1}, PerBaseCoverage{1, 0, 0, 0},
      PerBaseCoverage{0},       PerBaseCoverage{}};

  EXPECT_EQ(expected_coverage, actual_coverage);
}

TEST_F(PbCovRecorder_nestedDeletion,
       multiMappedReadSingleSearchState_correctRecordedPbCoverage) {
  // PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "CTTT"
  std::size_t read_size{4};

  PbCovRecorder{prg_info, multi_mapped_reads_2, read_size};
  auto actual_coverage =
      collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

  SitePbCoverage expected_coverage{
      PerBaseCoverage{},        PerBaseCoverage{0, 0},
      PerBaseCoverage{0, 0, 1}, PerBaseCoverage{0, 0, 0, 1},
      PerBaseCoverage{0},       PerBaseCoverage{}};

  EXPECT_EQ(expected_coverage, actual_coverage);
}
