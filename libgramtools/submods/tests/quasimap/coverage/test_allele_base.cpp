#include <cctype>

#include "gtest/gtest.h"


#include "src_common/generate_prg.hpp"
#include "quasimap/coverage/allele_base.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"


using namespace gram;
using namespace gram::coverage::per_base;

/*
PRG: GCT5C6AA6T6AG7T8C8CT
        i	BWT	SA	text_suffix
0	T	20	0
1	6	6	A A 6 T 6 A G 7 T 8 C 8 C T 0
2	6	11	A G 7 T 8 C 8 C T 0
3	A	7	A 6 T 6 A G 7 T 8 C 8 C T 0
4	8	18	C T 0
5	G	1	C T 5 C 6 A A 6 T 6 A G 7 T 8 C 8 C T 0
6	5	4	C 6 A A 6 T 6 A G 7 T 8 C 8 C T 0
7	8	16	C 8 C T 0
8	0	0	G C T 5 C 6 A A 6 T 6 A G 7 T 8 C 8 C T 0
9	A	12	G 7 T 8 C 8 C T 0
10	C	19	T 0
11	C	2	T 5 C 6 A A 6 T 6 A G 7 T 8 C 8 C T 0
12	6	9	T 6 A G 7 T 8 C 8 C T 0
13	7	14	T 8 C 8 C T 0
14	T	3	5 C 6 A A 6 T 6 A G 7 T 8 C 8 C T 0
15	C	5	6 A A 6 T 6 A G 7 T 8 C 8 C T 0
16	T	10	6 A G 7 T 8 C 8 C T 0
17	A	8	6 T 6 A G 7 T 8 C 8 C T 0
18	G	13	7 T 8 C 8 C T 0
19	C	17	8 C T 0
20	T	15	8 C 8 C T 0
*/

TEST(AlleleBaseCoverage, ReadCoversTwoSites_CorrectAlleleBaseCoverage) {
    auto prg_raw = encode_prg("gct5c6aa6t6ag7t8c8ct");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uint64_t read_length = 6;

    SearchState search_state = {
            SA_Interval{11, 11},
            VariantSitePath{
                    VariantLocus{7, 2},
                    VariantLocus{5, 2}
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{0}, {1, 1}, {0}},
            {{0}, {1}},
    };
    EXPECT_EQ(result, expected);
}

TEST(Site_Boundaries, Get_Start_Ends){
    auto prg_raw = encode_prg("gct5c6aa6t6ag7t8c8ct");
    auto prg_info = generate_prg_info(prg_raw);

    auto site_boundaries = gram::site_marker_prg_indexes(5,prg_info);
    EXPECT_EQ(site_boundaries.first, 3);
    EXPECT_EQ(site_boundaries.second, 10);
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

TEST(AlleleBaseCoverage, ShortReadStartingOutsideSiteCoversTwoSites_FinishesBeforeSecondAlleleEnd) {
    auto prg_raw = encode_prg("gct5c6g6t6ag7t8cc8ct");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uint64_t read_length = 6;

    SearchState search_state = {
            SA_Interval{4, 4},
            VariantSitePath{
                    VariantLocus{7, 2},
                    VariantLocus{5, 2},
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{0}, {1}, {0}},
            {{0}, {1, 0}},
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, ReadStartsWithinOneAlleleFinishesBeforeEndOfSecond_CorrectCoverage) {
    auto prg_raw = encode_prg("gct5c6g6t6ag7t8cc8ct");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval{12, 12},
            VariantSitePath{
                    VariantLocus{7, 2},
                    VariantLocus{5, 3},
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{0}, {0}, {1}},
            {{0}, {1, 0}},
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenTwoSites_CorrectInterSiteBaseCount) {
    auto prg_raw = encode_prg("gct5c6g6t6ag7t8cc8ct");
    auto prg_info = generate_prg_info(prg_raw);

    uint64_t first_site_marker = 5;
    uint64_t second_site_marker = 7;

    auto first_site_prg_start_end = site_marker_prg_indexes(first_site_marker, prg_info);
    auto first_site_prg_end = first_site_prg_start_end.second;

    auto second_site_prg_start_end = site_marker_prg_indexes(second_site_marker, prg_info);
    auto second_site_prg_start = second_site_prg_start_end.first;

    uint64_t result = second_site_prg_start - first_site_prg_end - 1;
    uint64_t expected = 2;
    EXPECT_EQ(result, expected);
}


TEST(SetSiteBaseCoverage, AlleleOffsetGreaterThanBasesToSet_CorrectBasesSet) {
    auto prg_raw = encode_prg("gct5c6agtaaatgcg6agt");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    SitesCoverageBoundaries sites_coverage_boundaries = {};
    VariantLocus path_element = {5, 2};
    uint64_t allele_coverage_offset = 6;
    uint64_t max_bases_to_set = 3;

    set_site_base_coverage(coverage,
                           sites_coverage_boundaries,
                           path_element,
                           allele_coverage_offset,
                           max_bases_to_set);

    const auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {
                    {0}, {0, 0, 0, 0, 0, 0, 1, 1, 1, 0}
            }
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: AC5GG6AGA6C
i	BWT	SA	text_suffix
0	C	11
1	0	0	A C 5 G G 6 A G A 6 C
2	6	6	A G A 6 C
3	G	8	A 6 C
4	6	10	C
5	A	1	C 5 G G 6 A G A 6 C
6	A	7	G A 6 C
7	5	3	G G 6 A G A 6 C
8	G	4	G 6 A G A 6 C
9	C	2	5 G G 6 A G A 6 C
10	G	5	6 A G A 6 C
11	A	9	6 C
*/

TEST(AlleleBaseCoverage, SaIntervalGreaterThanOne_CorrectCumulativeBaseCoverage) {
    auto prg_raw = encode_prg("ac5gg6aga6c");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval{7, 8},
            VariantSitePath{
                    VariantLocus{5, 1},
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{1, 1}, {0, 0, 0}}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, ReadStartsBeforeSiteCoversFirstAllele_CorrectBaseCoverage) {
    auto prg_raw = encode_prg("ac5gg6aga6c");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uint64_t read_length = 150;

    SearchState search_state = {
            SA_Interval{1, 1},
            VariantSitePath{
                    VariantLocus{5, 1}
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{1, 1}, {0, 0, 0}}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, ReadStartsWithinFirstAllele_OnlyLastAlleleBaseCovered) {
    auto prg_raw = encode_prg("ac5gg6aga6c");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uint64_t read_length = 150;

    SearchState search_state = {
            SA_Interval{8, 8},
            VariantSitePath{
                    VariantLocus{5, 1}
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{0, 1}, {0, 0, 0}}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, ReadStartsWithinSecondAllele_PartialAlleleBaseCoverage) {
    auto prg_raw = encode_prg("ac5gg6aga6c");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uint64_t read_length = 150;

    SearchState search_state = {
            SA_Interval{6, 6},
            VariantSitePath{
                    VariantLocus{5, 2}
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{0, 0}, {0, 1, 1}}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, ReadStartsOutsideSiteEndsBeforeAlleleEnd_PartialCoverageOfAllele) {
    auto prg_raw = encode_prg("ac5gg6aga6c");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval{1, 1},
            VariantSitePath{
                    VariantLocus{5, 2}
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{0, 0}, {1, 1, 0}}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenSiteStartingAtPrgStart_CorrectAlleleBaseCoverageStructure) {
    auto prg_raw = encode_prg("5gg6aga6c");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_base_structure(prg_info);
    SitesAlleleBaseCoverage expected = {
            AlleleCoverage{
                    BaseCoverage{0, 0},
                    BaseCoverage{0, 0, 0},
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenOneVariantSite_CorrectAlleleBaseCoverageStructure) {
    auto prg_raw = encode_prg("ct5gg6aga6c");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_base_structure(prg_info);
    SitesAlleleBaseCoverage expected = {
            AlleleCoverage{
                    BaseCoverage{0, 0},
                    BaseCoverage{0, 0, 0},
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenTwoVariantSites_CorrectAlleleBaseCoverageStructure) {
    auto prg_raw = encode_prg("ct5gg6aga6ccccc7a8ttt8");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_base_structure(prg_info);
    SitesAlleleBaseCoverage expected = {
            AlleleCoverage{
                    BaseCoverage{0, 0},
                    BaseCoverage{0, 0, 0},
            },
            AlleleCoverage{
                    BaseCoverage{0},
                    BaseCoverage{0, 0, 0},
            },
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenPopulatedAlleleBaseCoverage_CorrectJsonDump) {
    SitesAlleleBaseCoverage allele_base_coverage = {
            AlleleCoverage{
                    BaseCoverage{1, 12},
                    BaseCoverage{0, 3, 0},
            },
            AlleleCoverage{
                    BaseCoverage{0},
                    BaseCoverage{0, 19, 0},
            },
    };
    auto result = dump_allele_base_coverage(allele_base_coverage);
    std::string expected = "{\"allele_base_counts\":[[[1,12],[0,3,0]],[[0],[0,19,0]]]}";
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenSingleSiteAlleleBaseCoverage_CorrectJsonDump) {
    SitesAlleleBaseCoverage allele_base_coverage = {
            AlleleCoverage{
                    BaseCoverage{1, 12},
                    BaseCoverage{0, 3, 0},
            }
    };
    auto result = dump_allele_base_coverage(allele_base_coverage);
    std::string expected = "{\"allele_base_counts\":[[[1,12],[0,3,0]]]}";
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenEmptyAlleleBaseCoverage_CorrectJsonDump) {
    SitesAlleleBaseCoverage allele_base_coverage = {};
    auto result = dump_allele_base_coverage(allele_base_coverage);
    std::string expected = "{\"allele_base_counts\":[]}";
    EXPECT_EQ(result, expected);
}


TEST(AlleleStartOffsetIndex, GivenSecondAlleleBase_CorrectAlleleIndexOffset) {
    auto prg_raw = encode_prg("ct5gg6aaga5cc");
    auto prg_info = generate_prg_info(prg_raw);

    uint64_t within_allele_prg_index = 7;
    uint64_t result = allele_start_offset_index(within_allele_prg_index, prg_info);
    uint64_t expected = 1;

    EXPECT_EQ(expected, result);
}


TEST(AlleleStartOffsetIndex, GivenFirstAlleleBase_CorrectAlleleIndexOffset) {
    auto prg_raw = encode_prg("ct5gg6aaga5cc");
    auto prg_info = generate_prg_info(prg_raw);

    uint64_t within_allele_prg_index = 6;
    uint64_t result = allele_start_offset_index(within_allele_prg_index, prg_info);
    uint64_t expected = 0;

    EXPECT_EQ(expected, result);
}

TEST(DummyCovNode, BuildWithSizeSmallerThanEndCoord_ThrowsException){
   EXPECT_THROW(DummyCovNode(0, 5, 3), InconsistentCovNodeCoordinates);
}

TEST(DummyCovNode, BuildWithStartGreaterThanEnd_ThrowsException){
    EXPECT_THROW(DummyCovNode(2, 1, 3), InconsistentCovNodeCoordinates);
}

TEST(DummyCovNode, ExtendWithEndPosGreaterThanNodeSize_ThrowsException){
    auto d = DummyCovNode(1, 5, 6);
    EXPECT_THROW(d.extend_coordinates({0, 6}), InconsistentCovNodeCoordinates);
}

TEST(DummyCovNode, ExtendNoStartAndNoEnd_CorrectUnchangesCoordinates) {
    auto d = DummyCovNode(1, 5, 6);
    d.extend_coordinates({2, 5});
    node_coordinates expected_coords{1, 5};
    EXPECT_EQ(expected_coords, d.get_coordinates());
}

TEST(DummyCovNode, ExtendStartAndEnd_CorrectExtendedCoordinates){
    auto d = DummyCovNode(3, 3, 6);
    d.extend_coordinates({0, 5});
    node_coordinates expected_coords{0, 5};
    EXPECT_EQ(expected_coords, d.get_coordinates());
}

TEST(Traverser, StartOutOfSiteEndInSite_correctObjectState){
    auto prg_raw = encode_prg("CT5gg6AAGa5cc");
    auto prg_info = generate_prg_info(prg_raw);

    std::size_t read_size = 5;
    VariantSitePath traversed_path{
        VariantLocus{5, 2}
    };
    auto start_point = prg_info.coverage_graph.random_access[0];

    Traverser t{start_point, traversed_path, read_size};
    auto variant_node = t.next_Node().value();
    EXPECT_EQ(variant_node->get_site(), 5);
    EXPECT_EQ(variant_node->get_allele(), 2);

    std::pair<uint32_t,uint32_t> expected_coordinates{0, 2};
    EXPECT_EQ(expected_coordinates, t.get_node_coordinates());
    EXPECT_EQ(false, t.next_Node().has_value());
}

TEST(Traverser, StartAndEndInSite_CorrectNodeInterval){
    auto prg_raw = encode_prg("ct5g6aaAAAAAAaga5cc");
    auto prg_info = generate_prg_info(prg_raw);

    std::size_t read_size = 6;
    // Empty because the fact we are in VariantLocus{5, 2} is recorded in traversing_path container
    VariantSitePath traversed_path{};
    auto start_point = prg_info.coverage_graph.random_access[7];

    Traverser t{start_point, traversed_path, read_size};
    auto variant_node = t.next_Node().value();

    std::pair<uint32_t,uint32_t> expected_coordinates{2, 7};
    EXPECT_EQ(expected_coordinates, t.get_node_coordinates());
}

TEST(Traverser, StartInSiteAndTraverseToAnotherSite_CorrectObjectState){
    auto prg_raw = encode_prg("ct5g6aAA6CC7gc8ga8AAAAa8");
    auto prg_info = generate_prg_info(prg_raw);

    std::size_t read_size = 8;
    VariantSitePath traversed_path{
        VariantLocus{7, 3}
    };
    auto start_point = prg_info.coverage_graph.random_access[6];

    Traverser t{start_point, traversed_path, read_size};
    auto cur_Node = t.next_Node();
    auto variant_node = cur_Node;
    while (cur_Node.has_value()) {
        variant_node = cur_Node;
        cur_Node = t.next_Node();
    }

    std::pair<uint32_t,uint32_t> expected_coordinates{0, 3};
    EXPECT_EQ(expected_coordinates, t.get_node_coordinates());
    EXPECT_EQ(0, t.get_remaining_bases());
}

// Helper function to get all the loci that were traversed. Modifies the traversal in place
VariantSitePath collect_traversal(Traverser & t){
    VariantSitePath traversal;
    VariantLocus site_and_allele;
    auto cur_Node = t.next_Node();

    while(bool(cur_Node)){
        site_and_allele = {cur_Node.value()->get_site(), cur_Node.value()->get_allele()};
        traversal.push_back(site_and_allele);
        cur_Node = t.next_Node();
    }
    return traversal;
}

TEST(Traverser_Nested, StartOutOfSiteEndOutOfSite_CorrectChosenSitesAndEndState){
    std::string raw_prg = "A[ctt,G[AAA,a]T]CCccc";
    marker_vec v = prg_string_to_ints(raw_prg);
    auto prg_info = generate_prg_info(v);

    std::size_t read_size = 8;
    VariantSitePath traversed_path{
            VariantLocus{7, 1},
            VariantLocus{5, 2}
    };

    auto start_point = prg_info.coverage_graph.random_access[0];
    Traverser t{start_point, traversed_path, read_size};

    VariantSitePath expected_traversal{
        VariantLocus{5, 2},
        VariantLocus{7, 1},
        VariantLocus{5, 2} // After exiting site 7, we still have coverage to record on allele 2 of site 5 (base 'T')
    };

    VariantSitePath actual_traversal = collect_traversal(t);
    EXPECT_EQ(expected_traversal, actual_traversal);

    // Make sure we have consumed all bases of the read
    EXPECT_EQ(0, t.get_remaining_bases());
    // Make sure we are placed correctly in the last node
    std::pair<uint32_t, uint32_t> expected_last_node_coords{0, 1};
    EXPECT_EQ(expected_last_node_coords, t.get_node_coordinates());
}

TEST(Traverser_Nested, TraverseGraphWithLevel2Nesting_CorrectChosenSitesAndEndState){
    std::string raw_prg = "A[CT[GC[c,A]A,gc]T[C,a]Tt,t]c";
    marker_vec v = prg_string_to_ints(raw_prg);
    auto prg_info = generate_prg_info(v);

    std::size_t read_size = 10;
    VariantSitePath traversed_path{
            VariantLocus{11, 1},
            VariantLocus{9, 2},
            VariantLocus{7, 1},
            VariantLocus{5, 1}
    };
    auto start_point = prg_info.coverage_graph.random_access[0];
    Traverser t{start_point, traversed_path, read_size};

    VariantSitePath expected_traversal{
            VariantLocus{5, 1},
            VariantLocus{7, 1},
            VariantLocus{9, 2},
            VariantLocus{7, 1},
            VariantLocus{5, 1},
            VariantLocus{11, 1},
            VariantLocus{5, 1},
    };

    VariantSitePath actual_traversal = collect_traversal(t);
    EXPECT_EQ(expected_traversal, actual_traversal);

    EXPECT_EQ(0, t.get_remaining_bases());
    std::pair<uint32_t, uint32_t> expected_last_node_coords{0, 0};
    EXPECT_EQ(expected_last_node_coords, t.get_node_coordinates());
}

TEST(PbCovRecorder,ProcessNewCovNode_CorrectDummyCovNodeMade){
    PbCovRecorder pb_rec;
    covG_ptr cov_node = boost::make_shared<coverage_Node>(
            coverage_Node{"ACTG", 102, 5, 2});
    realCov_to_dummyCov expected_mapping{
            {cov_node, DummyCovNode(1, 3, 4)}
    };

    pb_rec.process_Node(cov_node, 1, 3);
    EXPECT_EQ(expected_mapping, pb_rec.get_cov_mapping());
}

TEST(PbCovRecorder, ProcessExistingCovNode_CorrectlyUpdatedDummyCovNode){
    covG_ptr cov_node = boost::make_shared<coverage_Node>(
            coverage_Node{"ACTGCC", 99, 5, 2});
    realCov_to_dummyCov existing_mapping  {
            {cov_node, DummyCovNode{1, 3, 6}}
    };
    PbCovRecorder pb_rec(existing_mapping);
    pb_rec.process_Node(cov_node, 2, 5);

    realCov_to_dummyCov expected_mapping{
            {cov_node, DummyCovNode(1, 5, 6)}
    };

    EXPECT_EQ(expected_mapping, pb_rec.get_cov_mapping());
}

using prg_positions = std::vector<std::size_t>;
using dummy_cov_nodes = std::vector<DummyCovNode>;
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
    SearchState simple_read_1{
            SA_Interval{5, 5},
            VariantSitePath{VariantLocus{7, 1}}
    };

    // Read: ATTTT
    SearchState simple_read_2{
            SA_Interval{1, 1},
            VariantSitePath{VariantLocus{5, 2}}
    };

    // Read: GCC. Two distinct occurrences compatible with same sites
    SearchStates multi_mapped_reads_1{
            SearchState{
                    SA_Interval{9, 9},
                    VariantSitePath{VariantLocus{7, 1}}
            },
            SearchState{
                    SA_Interval{8, 8},
                    VariantSitePath{}
            }
    };

    // Read: CTT. In a single Search State
    SearchStates multi_mapped_reads_2{
            SearchState{
                    SA_Interval{6, 7},
                    VariantSitePath{}
            }
    };
};


dummy_cov_nodes collect_dummy_cov_nodes(coverage_Graph const& cov_graph, prg_positions positions,
        realCov_to_dummyCov cov_mapping){
    dummy_cov_nodes result(positions.size());
    covG_ptr accessed_node;
    std::size_t index{0};

    for (auto& pos : positions){
       accessed_node = cov_graph.random_access[pos].node;
       if (cov_mapping.find(accessed_node) == cov_mapping.end()) result[index] = DummyCovNode{};
       else result[index] = cov_mapping.at(accessed_node);
        index++;
    }
    return result;
}


AlleleCoverage collect_coverage(coverage_Graph const& cov_graph, prg_positions positions){
   AlleleCoverage result(positions.size());
    covG_ptr accessed_node;
    std::size_t index{0};

    for (auto& pos : positions) {
        accessed_node = cov_graph.random_access[pos].node;
        result[index] = accessed_node->get_coverage();
        index++;
    }
   return result;
}

TEST_F(PbCovRecorder_nestedDeletion, simpleRead1Mapped_correctDummyCovNodes){
   //PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "CGCCTT"
   std::size_t read_size{6};
   PbCovRecorder recorder(prg_info, read_size);
   recorder.process_SearchState(simple_read_1);
   auto actual_dummies = collect_dummy_cov_nodes(prg_info.coverage_graph, all_sequence_node_positions,
           recorder.get_cov_mapping());

   dummy_cov_nodes expected_dummies{
           DummyCovNode{}, DummyCovNode{1, 1, 2},
           DummyCovNode{0, 2, 3}, DummyCovNode{},
           DummyCovNode{}, DummyCovNode{}
   };
   EXPECT_EQ(expected_dummies, actual_dummies);
}

TEST_F(PbCovRecorder_nestedDeletion, simpleRead1Mapped_correctRecordePbCoverage){
    //PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "CGCCTT"
    SearchStates mapping{simple_read_1};
    std::size_t read_size{6};
    PbCovRecorder recorder(prg_info, mapping, read_size);
    auto actual_coverage = collect_coverage(prg_info.coverage_graph, all_sequence_node_positions);

    AlleleCoverage expected_coverage{
            BaseCoverage{}, BaseCoverage{0, 1},
            BaseCoverage{1, 1, 1}, BaseCoverage{0, 0, 0, 0},
            BaseCoverage{0}, BaseCoverage{}
    };
    EXPECT_EQ(expected_coverage, actual_coverage);
}

TEST_F(PbCovRecorder_nestedDeletion, simpleRead2Mapped_correctDummyCovNodes){
    //PRG: "AT[GC[GCC,CCGC],T]TTTT"; Read: "ATTTT"
    SearchStates mapping = SearchStates{simple_read_2};
    std::size_t read_size{5};
    PbCovRecorder recorder(prg_info, read_size);
    recorder.process_SearchState(simple_read_2);
    auto actual_dummies = collect_dummy_cov_nodes(prg_info.coverage_graph, all_sequence_node_positions,
                                                  recorder.get_cov_mapping());

    dummy_cov_nodes expected_dummies{
            DummyCovNode{}, DummyCovNode{},
            DummyCovNode{}, DummyCovNode{},
            DummyCovNode{0, 0, 1}, DummyCovNode{}
    };
    EXPECT_EQ(expected_dummies, actual_dummies);
}
