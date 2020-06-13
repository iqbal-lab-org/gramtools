#include "gtest/gtest.h"

#include "genotype/quasimap/coverage/coverage_common.hpp"

#include "../../../test_resources/mocks.hpp"
#include "submod_resources.hpp"

using namespace gram;

std::set<SitePath> get_site_path_only(uniqueSitePaths const& map) {
  std::set<SitePath> site_path;
  for (auto const& e : map) {
    site_path.insert(e.first);
  }
  return site_path;
}

TEST(RandomInclusiveInt, RandomCall_MinBoundaryReturned) {
  uint32_t random_seed = 48;
  RandomInclusiveInt r{random_seed};
  uint32_t result = r.generate(1, 10);
  uint32_t expected = 5;
  EXPECT_EQ(result, expected);
}

TEST(RandomInclusiveInt, RandomCall_MaxBoundaryReturned) {
  uint32_t random_seed = 56;
  RandomInclusiveInt r{random_seed};
  uint32_t result = r.generate(1, 10);
  uint32_t expected = 4;
  EXPECT_EQ(result, expected);
}

TEST(CountNonvariantSearchStates, OnePathOneNonPath_CountOne) {
  SearchStates search_states = {
      SearchState{SA_Interval{},
                  VariantSitePath{
                      VariantLocus{5, FIRST_ALLELE},
                      VariantLocus{7, FIRST_ALLELE + 1},
                  }},
      SearchState{SA_Interval{}, VariantSitePath{}}

  };
  MappingInstanceSelector s;
  auto result = s.count_nonvar_search_states(search_states);
  uint64_t expected = 1;
  EXPECT_EQ(result, expected);
}

TEST(LocusFinder_logic, SameSiteMoreThanOnceInSearchState_ThrowsError) {
  SearchState search_state = SearchState{
      SA_Interval{}, VariantSitePath{VariantLocus{5, FIRST_ALLELE + 1}},
      VariantSitePath{VariantLocus{5, ALLELE_UNKNOWN}}};
  LocusFinder l;
  EXPECT_THROW(l.check_site_uniqueness(search_state), std::logic_error);
}

TEST(SameLevel0SitesDifferentOrder, SingleEntryInMap) {
  level0_Sites s1{Marker{5}, Marker{7}, Marker{9}, Marker{11}};
  level0_Sites s2{Marker{11}, Marker{9}, Marker{7}, Marker{5}};

  uniqueSitePaths unique_map{};
  unique_map.insert({s1, traversal_info{}});
  unique_map.insert({s2, traversal_info{}});

  EXPECT_EQ(1, unique_map.size());
}

TEST(GetUniquePathSites, TwoDifferentPaths_CorrectPaths) {
  SearchStates search_states = {
      SearchState{SA_Interval{},
                  VariantSitePath{
                      VariantLocus{5, FIRST_ALLELE},
                      VariantLocus{7, FIRST_ALLELE + 1},
                  }},
      SearchState{SA_Interval{}, VariantSitePath{
                                     VariantLocus{9, FIRST_ALLELE + 2},
                                     VariantLocus{11, FIRST_ALLELE + 4},
                                 }}};
  PRG_Info prg_info;
  MappingInstanceSelector m{&prg_info};
  m.process_searchstates(search_states);
  auto result = get_site_path_only(m.usps);
  std::set<SitePath> expected = {SitePath{5, 7}, SitePath{9, 11}};
  EXPECT_EQ(result, expected);

  // Check SearchState dispatch
  auto ss_result1 = m.usps[SitePath{5, 7}].first.front();
  EXPECT_EQ(ss_result1, search_states.front());

  auto ss_result2 = m.usps[SitePath{9, 11}].first.front();
  EXPECT_EQ(ss_result2, search_states.back());
}

TEST(GetUniquePathSites,
     TwoIdenticalPathsOneEmptyPath_SingleNonEmptyPathInSet) {
  SearchStates search_states = {
      SearchState{SA_Interval{},
                  VariantSitePath{
                      VariantLocus{9, FIRST_ALLELE + 2},
                      VariantLocus{11, FIRST_ALLELE + 4},
                  }},
      SearchState{SA_Interval{},
                  VariantSitePath{
                      VariantLocus{9, FIRST_ALLELE + 2},
                      VariantLocus{11, FIRST_ALLELE + 4},
                  }},
      SearchState{SA_Interval{}, VariantSitePath{}}

  };
  PRG_Info prg_info;
  MappingInstanceSelector m{&prg_info};
  m.process_searchstates(search_states);
  auto result = get_site_path_only(m.usps);
  std::set<SitePath> expected = {SitePath{9, 11}};
  EXPECT_EQ(result, expected);
}

class LocusFinder_minimal : public ::testing::Test {
 protected:
  void SetUp() {
    parental_map p{{9, VariantLocus{7, FIRST_ALLELE}},
                   {7, VariantLocus{5, FIRST_ALLELE + 2}}};
    coverage_Graph c;
    c.par_map = p;
    prg_info.coverage_graph = c;
  }
  LocusFinder l{};
  PRG_Info prg_info;
};

TEST_F(LocusFinder_minimal, assignNestedLocus_correctDispatching) {
  // First addition
  VariantLocus test{9, FIRST_ALLELE + 2};
  l.assign_nested_locus(test, &prg_info);
  SitePath expected_base_sites{5};
  EXPECT_EQ(l.base_sites, expected_base_sites);

  SitePath expected_used_sites{5, 7, 9};
  EXPECT_EQ(l.used_sites, expected_used_sites);

  uniqueLoci expected_unique_loci{VariantLocus{5, FIRST_ALLELE + 2},
                                  VariantLocus{7, FIRST_ALLELE},
                                  VariantLocus{9, FIRST_ALLELE + 2}};
  EXPECT_EQ(l.unique_loci, expected_unique_loci);

  // Second addition: nothing should change
  VariantLocus test2{9, 2};
  l.assign_nested_locus(test2, &prg_info);
  EXPECT_EQ(l.base_sites, expected_base_sites);
  EXPECT_EQ(l.used_sites, expected_used_sites);
  EXPECT_EQ(l.unique_loci, expected_unique_loci);
}

TEST_F(LocusFinder_minimal, assignTraversedLoci_correctDispatching) {
  SearchState test{SA_Interval{2, 2},
                   VariantSitePath{VariantLocus{11, FIRST_ALLELE},
                                   VariantLocus{9, FIRST_ALLELE + 2}}};

  l.assign_traversed_loci(test, &prg_info);
  SitePath expected_base_sites{5, 11};
  EXPECT_EQ(l.base_sites, expected_base_sites);

  uniqueLoci expected_unique_loci{
      VariantLocus{5, FIRST_ALLELE + 2}, VariantLocus{7, FIRST_ALLELE},
      VariantLocus{9, FIRST_ALLELE + 2}, VariantLocus{11, FIRST_ALLELE}};
  EXPECT_EQ(l.unique_loci, expected_unique_loci);
}

/*
PRG: A[[G[AC,TC],A]C,T]T
i	BWT	SA	text_suffix
0	T	19	0
1	9	5	A C 10 T C 10 8 A 8 C 6 T 6 T 0
2	0	0	A 5 7 G 9 A C 10 T C 10 8 A 8 C 6 T 6 T 0
3	8	12	A 8 C 6 T 6 T 0
4	8	14	C 6 T 6 T 0
5	A	6	C 10 T C 10 8 A 8 C 6 T 6 T 0
6	T	9	C 10 8 A 8 C 6 T 6 T 0
7	7	3	G 9 A C 10 T C 10 8 A 8 C 6 T 6 T 0
8	6	18	T 0
9	10	8	T C 10 8 A 8 C 6 T 6 T 0
10	6	16	T 6 T 0
11	A	1	5 7 G 9 A C 10 T C 10 8 A 8 C 6 T 6 T 0
12	T	17	6 T 0
13	C	15	6 T 6 T 0
14	5	2	7 G 9 A C 10 T C 10 8 A 8 C 6 T 6 T 0
15	10	11	8 A 8 C 6 T 6 T 0
16	A	13	8 C 6 T 6 T 0
17	G	4	9 A C 10 T C 10 8 A 8 C 6 T 6 T 0
18	C	7	10 T C 10 8 A 8 C 6 T 6 T 0
19	C	10	10 8 A 8 C 6 T 6 T 0 */

class LocusFinder_full : public ::testing::Test {
  /**
   * Here we make a full fm index and coverage graph
   * Disclaimer: the tests are strongly coupled to, and thus require correctness
   * of: i) Coverage graph (parent_map; random_access to nodes) ii) FM Index
   * construction We could [TODO] decouple and write/mock those ourselves.
   */
  //
 protected:
  void SetUp() {
    std::string raw_prg = "A[[G[AC,TC],A]C,T]T";
    marker_vec v = prg_string_to_ints(raw_prg);
    prg_info = generate_prg_info(v);
  }
  LocusFinder l;
  PRG_Info prg_info;
};

TEST_F(LocusFinder_full,
       assignTraversingLociWithAllUnknownLoci_correctDispatching) {
  // Pretense is we've mapped the read "CCT"
  SearchState test{SA_Interval{5, 6}, VariantSitePath{},
                   VariantSitePath{VariantLocus{5, ALLELE_UNKNOWN},
                                   VariantLocus{7, ALLELE_UNKNOWN},
                                   VariantLocus{9, ALLELE_UNKNOWN}}};
  l.assign_traversing_loci(test, &prg_info);

  SitePath expected_base_sites{5};
  EXPECT_EQ(l.base_sites, expected_base_sites);

  uniqueLoci expected_unique_loci{
      VariantLocus{5, FIRST_ALLELE}, VariantLocus{7, FIRST_ALLELE},
      VariantLocus{9, FIRST_ALLELE}, VariantLocus{9, FIRST_ALLELE + 1}};
  EXPECT_EQ(l.unique_loci, expected_unique_loci);
}

TEST_F(LocusFinder_full,
       assignTraversedLociWithOneTraversedLocus_correctDispatching) {
  // Pretense is we've mapped the read "GACC"
  SearchState test{SA_Interval{7, 7},
                   VariantSitePath{VariantLocus{9, FIRST_ALLELE}},
                   VariantSitePath{
                       VariantLocus{7, ALLELE_UNKNOWN},
                   }};

  l.assign_traversing_loci(test, &prg_info);
  SitePath expected_base_sites{5};
  EXPECT_EQ(l.base_sites, expected_base_sites);

  uniqueLoci expected_unique_loci{
      VariantLocus{5, FIRST_ALLELE},
      VariantLocus{7, FIRST_ALLELE},
  };
  EXPECT_EQ(l.unique_loci, expected_unique_loci);
}

TEST_F(LocusFinder_full,
       constructLocusFinder_assignAllLociForSearchState_correctDispatching) {
  // Pretense is we've mapped the read "GACC"
  SearchState test{SA_Interval{7, 7},
                   VariantSitePath{VariantLocus{9, FIRST_ALLELE}},
                   VariantSitePath{
                       VariantLocus{7, ALLELE_UNKNOWN},
                   }};
  LocusFinder l2{test, &prg_info};

  SitePath expected_base_sites{5};
  EXPECT_EQ(l2.base_sites, expected_base_sites);

  uniqueLoci expected_unique_loci{
      VariantLocus{5, FIRST_ALLELE},
      VariantLocus{7, FIRST_ALLELE},
      VariantLocus{9, FIRST_ALLELE},
  };
  EXPECT_EQ(l2.unique_loci, expected_unique_loci);
}

class MappingInstanceSelector_addSearchStates : public ::testing::Test {
 protected:
  // In this example we pretend we have mapped "TAA" to the graph.
  // Note: the allele encapsulated mapping handling has separated a single
  // SearchState into three. The SA_Intervals are dummies.
  void SetUp() {
    std::string prg_raw{"[CG[TAA,T],TAA]TA[TAA,ATA]"};
    coverage_Graph c;
    parental_map par_map{{7, VariantLocus{5, FIRST_ALLELE}}};
    c.par_map = par_map;
    prg_info.coverage_graph = c;
  };
  PRG_Info prg_info;
  SearchState s1{SA_Interval{1, 1},
                 VariantSitePath{VariantLocus{7, FIRST_ALLELE}}};

  SearchState s2{SA_Interval{1, 1},
                 VariantSitePath{VariantLocus{5, FIRST_ALLELE + 1}}};
  SearchState s3{SA_Interval{1, 1},
                 VariantSitePath{VariantLocus{9, FIRST_ALLELE}}};
  MappingInstanceSelector selector{&prg_info};
};

TEST_F(MappingInstanceSelector_addSearchStates,
       addOneSearchState_correctlyRegistered) {
  selector.add_searchstate(s1);
  traversal_info expected_info{
      SearchStates{s1},
      uniqueLoci{VariantLocus{5, FIRST_ALLELE}, VariantLocus{7, FIRST_ALLELE}}};
  uniqueSitePaths expected_map{{SitePath{5}, expected_info}};

  EXPECT_EQ(selector.usps, expected_map);
}

TEST_F(MappingInstanceSelector_addSearchStates,
       addAllSearchStates_correctlyRegistered) {
  SearchStates all_ss{s1, s2, s3};
  selector.process_searchstates(all_ss);

  traversal_info expected_i1{
      SearchStates{s1, s2},
      uniqueLoci{VariantLocus{5, FIRST_ALLELE}, VariantLocus{7, FIRST_ALLELE},
                 VariantLocus{5, FIRST_ALLELE + 1}}};

  traversal_info expected_i2{SearchStates{s3},
                             uniqueLoci{VariantLocus{9, FIRST_ALLELE}}};

  uniqueSitePaths expected_map{{SitePath{5}, expected_i1},
                               {SitePath{9}, expected_i2}};
  EXPECT_EQ(selector.usps, expected_map);
}

class MappingInstanceSelector_select : public ::testing::Test {
  /*
   * There are four `SearchState`s: two go through two alleles of the same site
   * (7), and two do not cross any variant site in the PRG, and are held
   * together (SA Interval of size 2)
   *
   * The logic for random selection is choosing between 1 and 3, where 1 and 2
   * correspond to the invariants.
   */
 protected:
  PRG_Info prg_info;
  SearchStates ss{
      SearchState{SA_Interval{1, 1},
                  VariantSitePath{VariantLocus{7, FIRST_ALLELE}}},
      SearchState{SA_Interval{6, 6},
                  VariantSitePath{VariantLocus{7, FIRST_ALLELE + 1}}},
      SearchState{SA_Interval{2, 3}, VariantSitePath{}, VariantSitePath{}}};
};

TEST_F(MappingInstanceSelector_select,
       SelectInvariantAndThenVariant_correctIndices) {
  using namespace ::testing;
  MockRandomGenerator r;
  EXPECT_CALL(r, generate(1, 3))
      .Times(Exactly(2))
      .WillOnce(Return(1))
      .WillOnce(Return(3));

  MappingInstanceSelector m{&prg_info, &r};
  m.set_searchstates(ss);  // Copies them to m
  m.process_searchstates(ss);
  EXPECT_EQ(m.usps.size(), 1);  // Expect one unique site recorded: 7

  int32_t selected_index;
  selected_index = m.random_select_entry();
  EXPECT_EQ(selected_index,
            -1);  // We chose an invariant index, which returns -1

  selected_index = m.random_select_entry();
  EXPECT_EQ(
      selected_index,
      0);  // We chose the `SearchState` overlapping a variant site, at index 0.
}

TEST_F(MappingInstanceSelector_select, selectnonvariant_emptyMappingSelector) {
  // Select the SearchState in invariant region of PRG
  // The SA Interval is size 2 so the first two choices map to invariant mapping
  // instances
  using namespace ::testing;
  MockRandomGenerator r;
  EXPECT_CALL(r, generate(1, 3)).Times(Exactly(1)).WillOnce(Return(1));

  MappingInstanceSelector m{ss, &prg_info, &r};
  auto selection = m.get_selection();

  EXPECT_EQ(selection.navigational_search_states.size(), 0);
  EXPECT_EQ(selection.equivalence_class_loci.size(), 0);
}

TEST_F(MappingInstanceSelector_select, selectvariant_nonemptyMappingSelector) {
  using namespace ::testing;
  MockRandomGenerator r;
  EXPECT_CALL(r, generate(1, 3)).Times(Exactly(1)).WillOnce(Return(3));

  MappingInstanceSelector m{ss, &prg_info, &r};
  auto selection = m.get_selection();
  EXPECT_EQ(selection.navigational_search_states.size(), 2);
  uniqueLoci expected_loci{{VariantLocus{7, FIRST_ALLELE}},
                           {VariantLocus{7, FIRST_ALLELE + 1}}};
  EXPECT_EQ(selection.equivalence_class_loci, expected_loci);
}
