#include <cctype>
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "src_common/generate_prg.hpp"
#include "quasimap/coverage/common.hpp"

using namespace gram;

std::set<SitePath> get_site_path_only(uniqueSitePaths const& map){
    std::set<SitePath> site_path;
    for (auto const& e : map){
        site_path.insert(e.first);
    }
    return site_path;
}

/*
PRG: AA5T6CAGTAGCAGT6TA
i	BWT	SA	text_suffix
0	A	18
1	T	17	A
2	0	0	A A 5 T 6 C A G T A G C A G T 6 T A
3	T	9	A G C A G T 6 T A
4	C	6	A G T A G C A G T 6 T A
5	C	12	A G T 6 T A
6	A	1	A 5 T 6 C A G T A G C A G T 6 T A
7	6	5	C A G T A G C A G T 6 T A
8	G	11	C A G T 6 T A
9	A	10	G C A G T 6 T A
10	A	7	G T A G C A G T 6 T A
11	A	13	G T 6 T A
12	6	16	T A
13	G	8	T A G C A G T 6 T A
14	5	3	T 6 C A G T A G C A G T 6 T A
15	G	14	T 6 T A
16	A	2	5 T 6 C A G T A G C A G T 6 T A
17	T	4	6 C A G T A G C A G T 6 T A
18	T	15	6 T A
*/

TEST(CheckAlleleEncapsulated, TwoAlleleEncapsulatedMappings_True) {
    auto prg_raw = encode_prg("aa5t6cagtagcagt6ta");
    auto prg_info = generate_prg_info(prg_raw);

    // read: cagt
    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {7, 8},
            VariantSitePath {
                    VariantLocus {5, 2},
            },
            VariantSitePath {},
            SearchVariantSiteState::within_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_TRUE(result);
}


TEST(CheckAlleleEncapsulated, OneAlleleEncapsulatedMapping_True) {
    auto prg_raw = encode_prg("aa5t6cagtagcagt6ta");
    auto prg_info = generate_prg_info(prg_raw);

    // read: cagt
    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {7, 7},
            VariantSitePath {
                    VariantLocus {5, 2},
            },
            VariantSitePath {},
            SearchVariantSiteState::within_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_TRUE(result);
}


TEST(CheckAlleleEncapsulated, ReadOutsideOfSite_False) {
    auto prg_raw = encode_prg("aa5t6cagtagcagt6ta");
    auto prg_info = generate_prg_info(prg_raw);

    // read: aa
    uint64_t read_length = 2;

    SearchState search_state = {
            SA_Interval {2, 2},
            VariantSitePath {},
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_FALSE(result);
}


TEST(CheckAlleleEncapsulated, MappingExtendsOneBaseRightOustideOfSite_False) {
    auto prg_raw = encode_prg("aa5t6cagtagcAgt6ta");
    auto prg_info = generate_prg_info(prg_raw);

    // read: agtt
    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {5, 5},
            VariantSitePath {
                    VariantLocus {5, 2},
            },
            VariantSitePath {},
            SearchVariantSiteState::within_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_FALSE(result);
}


TEST(CheckAlleleEncapsulated, MappingExtendsOneBaseLeftOustideOfSite_False) {
    auto prg_raw = encode_prg("aa5t6cagtagcagt6ta");
    auto prg_info = generate_prg_info(prg_raw);

    // read: aca
    uint64_t read_length = 3;

    SearchState search_state = {
            SA_Interval {6, 6},
            VariantSitePath {
                    VariantLocus {5, 2},
            },
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_FALSE(result);
}


TEST(RandomInclusiveInt, RandomCall_MinBoundaryReturned) {
    uint32_t random_seed = 48;
    RandomInclusiveInt r{random_seed};
    uint32_t result = r.generate(1, 10);
    uint32_t expected = 1;
    EXPECT_EQ(result, expected);
}


TEST(RandomInclusiveInt, RandomCall_MaxBoundaryReturned) {
    uint32_t random_seed = 56;
    RandomInclusiveInt r{random_seed};
    uint32_t result = r.generate(1, 10);
    uint32_t expected = 10;
    EXPECT_EQ(result, expected);
}


TEST(CountNonvariantSearchStates, OnePathOneNonPath_CountOne) {
    SearchStates search_states = {
            SearchState {
                    SA_Interval {},
                    VariantSitePath {
                            VariantLocus {5, 1},
                            VariantLocus {7, 2},
                    }
            },
            SearchState {
                    SA_Interval {},
                    VariantSitePath {}
            }

    };
    MappingInstanceSelector s;
    auto result = s.count_nonvar_search_states(search_states);
    uint64_t expected = 1;
    EXPECT_EQ(result, expected);
}

TEST(GetSitePath, SameSiteMoreThanOnceInSearchState_ThrowsError){
    SearchState search_state = SearchState{
        SA_Interval{},
        VariantSitePath{
            VariantLocus{5, 2}
        },
        VariantSitePath{
            VariantLocus{5, ALLELE_UNKNOWN}
        }
    };

    EXPECT_THROW(get_path_sites(search_state), std::logic_error);
}

TEST(GetUniquePathSites, TwoDifferentPaths_CorrectPaths) {
    SearchStates search_states = {
            SearchState {
                    SA_Interval {},
                    VariantSitePath {
                            VariantLocus {5, 1},
                            VariantLocus {7, 2},
                    }
            },
            SearchState {
                    SA_Interval {},
                    VariantSitePath {
                            VariantLocus {9, 3},
                            VariantLocus {11, 5},
                    }
            }

    };
    auto result_map = get_unique_site_paths(search_states);
    auto result = get_site_path_only(result_map);
    std::set<SitePath> expected = {
            SitePath {5, 7},
            SitePath {9, 11}
    };
    EXPECT_EQ(result, expected);
}


TEST(GetUniquePathSites, TwoIdenticalPaths_SinglePathInSet) {
    SearchStates search_states = {
            SearchState {
                    SA_Interval {},
                    VariantSitePath {
                            VariantLocus {9, 3},
                            VariantLocus {11, 5},
                    }
            },
            SearchState {
                    SA_Interval {},
                    VariantSitePath {
                            VariantLocus {9, 3},
                            VariantLocus {11, 5},
                    }
            }

    };
    auto result_map = get_unique_site_paths(search_states);
    auto result = get_site_path_only(result_map);
    std::set<SitePath> expected = {
            SitePath {9, 11}
    };
    EXPECT_EQ(result, expected);
}


TEST(GetUniquePathSites, TwoIdenticalPathsOneEmptyPath_SingleNonEmptyPathInSet) {
    SearchStates search_states = {
            SearchState {
                    SA_Interval {},
                    VariantSitePath {
                            VariantLocus {9, 3},
                            VariantLocus {11, 5},
                    }
            },
            SearchState {
                    SA_Interval {},
                    VariantSitePath {
                            VariantLocus {9, 3},
                            VariantLocus {11, 5},
                    }
            },
            SearchState {
                    SA_Interval {},
                    VariantSitePath {}
            }

    };
    auto result_map = get_unique_site_paths(search_states);
    auto result = get_site_path_only(result_map);
    std::set<SitePath> expected = {
            SitePath {9, 11}
    };
    EXPECT_EQ(result, expected);
}


TEST(GetUniquePathSites, TwoSearchStatesSameSitePaths_CorrectUniquePathMap) {
    SearchStates search_states = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantLocus {5, 1},
                            VariantLocus {7, 2},
                    }
            },
            SearchState {
                    SA_Interval {3, 4},
                    VariantSitePath {
                            VariantLocus {5, 3},
                            VariantLocus {7, 2},
                    }
            }
    };

    uniqueSitePaths expected;
    expected.insert(std::make_pair(SitePath{5, 7}, search_states));
    auto result = get_unique_site_paths(search_states);
    EXPECT_EQ(result, expected);
}


TEST(GetUniquePathSites, SearchStatesWithSameAndDifferentSitePaths_CorrectUniquePathMap) {
    SearchStates same_search_states = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantLocus {5, 1},
                            VariantLocus {7, 2},
                    }
            },

            SearchState {
                    SA_Interval {5, 12},
                    VariantSitePath {
                            VariantLocus {5, 3},
                            VariantLocus {7, 5},
                    }
            }
    };

    SearchStates different_search_state = {
            SearchState {
                SA_Interval {3, 4},
                VariantSitePath {
                    VariantLocus {9, 3},
                    VariantLocus {11, 5},
                    }
            }
    };
    uniqueSitePaths expected;
    expected.insert(std::make_pair(SitePath{5, 7}, same_search_states));
    expected.insert(std::make_pair(SitePath{9, 11}, different_search_state));

    SearchStates all_search_states;
    all_search_states.insert(all_search_states.end(), same_search_states.begin(), same_search_states.end());
    all_search_states.emplace_back(different_search_state.front());
    auto result = get_unique_site_paths(all_search_states);

    EXPECT_EQ(result, expected);
}

class LocusFinder_minimal : public ::testing::Test{
protected:
    void SetUp(){
        parental_map p{
                {9, VariantLocus{7, 1}} ,
                {7, VariantLocus{5, 3}}
        };
        coverage_Graph c;
        c.par_map = p;
        prg_info.coverage_graph = c;
    }
    LocusFinder l{};
    PRG_Info prg_info;
};

TEST_F(LocusFinder_minimal, assignNestedLocus_correctDispatching){

    // First addition
    VariantLocus test{9, 3};
    l.assign_nested_locus(test, &prg_info);
    SitePath expected_base_sites{5};
    EXPECT_EQ(l.base_sites, expected_base_sites);

    SitePath expected_used_sites{5, 7, 9};
    EXPECT_EQ(l.used_sites, expected_used_sites);

    uniqueLoci expected_unique_loci{
        VariantLocus{5,3},
        VariantLocus{7,1},
        VariantLocus{9,3}
    };
    EXPECT_EQ(l.unique_loci, expected_unique_loci);

    //Second addition: nothing should change
    VariantLocus test2{9, 2};
    l.assign_nested_locus(test2, &prg_info);
    EXPECT_EQ(l.base_sites, expected_base_sites);
    EXPECT_EQ(l.used_sites, expected_used_sites);
    EXPECT_EQ(l.unique_loci, expected_unique_loci);
}

TEST_F(LocusFinder_minimal, assignTraversedLoci_correctDispatching){
    SearchState test{
        SA_Interval{2,2},
        VariantSitePath{
            VariantLocus{11, 1},
            VariantLocus{9, 3}
        }
    };

    l.assign_traversed_loci(test, &prg_info);
    SitePath expected_base_sites{5, 11};
    EXPECT_EQ(l.base_sites, expected_base_sites);

    uniqueLoci expected_unique_loci{
            VariantLocus{5,3},
            VariantLocus{7,1},
            VariantLocus{9,3},
            VariantLocus{11, 1}
    };
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

class LocusFinder_full : public ::testing::Test{
/**
 * Here we make a full fm index and coverage graph
 * Disclaimer: the tests are strongly coupled to, and thus require correctness of:
 *      i) Coverage graph (parent_map; random_access to nodes)
 *      ii) FM Index construction
 * We could [TODO] decouple and write/mock those ourselves.
 */
//
protected:
    void SetUp(){
       std::string raw_prg = "A[[G[AC,TC],A]C,T]T";
        marker_vec v = prg_string_to_ints(raw_prg);
        prg_info = generate_prg_info(v);
    }
    LocusFinder l;
    PRG_Info prg_info;
};

TEST_F(LocusFinder_full, assignTraversingLociWithAllUnknownLoci_correctDispatching){
    // Pretense is we've mapped the read "CCT"
    SearchState test{
            SA_Interval{5,6},
            VariantSitePath{},
            VariantSitePath{
                    VariantLocus{5, ALLELE_UNKNOWN},
                    VariantLocus{7,  ALLELE_UNKNOWN},
                    VariantLocus{9,  ALLELE_UNKNOWN}
            }
    };
    l.assign_traversing_loci(test, &prg_info);

    SitePath expected_base_sites{5};
    EXPECT_EQ(l.base_sites, expected_base_sites);

    uniqueLoci expected_unique_loci{
            VariantLocus{5,1},
            VariantLocus{7,1},
            VariantLocus{9,1},
            VariantLocus{9,2}
    };
    EXPECT_EQ(l.unique_loci, expected_unique_loci);
}

TEST_F(LocusFinder_full, assignTraversedLociWithOneTraversedLocus_correctDispatching) {
    // Pretense is we've mapped the read "GACC"
    SearchState test{
            SA_Interval{7,7},
            VariantSitePath{VariantLocus{9, 1}},
            VariantSitePath{
                    VariantLocus{7,  ALLELE_UNKNOWN},
            }
    };

    l.assign_traversing_loci(test, &prg_info);
    SitePath expected_base_sites{5};
    EXPECT_EQ(l.base_sites, expected_base_sites);

    uniqueLoci expected_unique_loci{
            VariantLocus{5,1},
            VariantLocus{7,1},
    };
    EXPECT_EQ(l.unique_loci, expected_unique_loci);
}

TEST_F(LocusFinder_full, constructLocusFinder_assignAllLociForSearchState_correctDispatching) {
    // Pretense is we've mapped the read "GACC"
    SearchState test{
            SA_Interval{7, 7},
            VariantSitePath{VariantLocus{9, 1}},
            VariantSitePath{
                    VariantLocus{7, ALLELE_UNKNOWN},
            }
    };
    LocusFinder l2{test, &prg_info};

    SitePath expected_base_sites{5};
    EXPECT_EQ(l2.base_sites, expected_base_sites);

    uniqueLoci expected_unique_loci{
            VariantLocus{5,1},
            VariantLocus{7,1},
            VariantLocus{9,1},
    };
    EXPECT_EQ(l2.unique_loci, expected_unique_loci);
}

class MappingInstanceSelector_addSearchStates : public ::testing::Test{
protected:
    // In this example we pretend we have mapped "TAA" to the graph.
    // Note: the allele encapsulated mapping handling has separated a single SearchState into three.
    // The SA_Intervals are dummies.
    void SetUp(){
        std::string prg_raw{"[CG[TAA,T],TAA]TA[TAA,ATA]"};
        coverage_Graph c;
        parental_map par_map{
                {7 , VariantLocus{5, 1}}
        };
        c.par_map = par_map;
        prg_info.coverage_graph = c;
    };
    PRG_Info prg_info;
    SearchState s1{
            SA_Interval{1, 1},
            VariantSitePath{VariantLocus{7, 1}}
    };

    SearchState s2{
            SA_Interval{1, 1},
            VariantSitePath{VariantLocus{5, 2}}
    };
    SearchState s3{
            SA_Interval{1, 1},
            VariantSitePath{VariantLocus{9, 1}}
    };
    MappingInstanceSelector selector{&prg_info};
};

TEST_F(MappingInstanceSelector_addSearchStates, addOneSearchState_correctlyRegistered){
    selector.add_searchstate(s1);
    traversal_info expected_info{
            SearchStates{s1},
            uniqueLoci{ VariantLocus{ 5, 1}, VariantLocus{7, 1} }
    };
    new_uniqueSitePaths expected_map{
            {SitePath{5} ,  expected_info}
    };

    EXPECT_EQ(selector.usps, expected_map);
}


TEST_F(MappingInstanceSelector_addSearchStates, addAllSearchStates_correctlyRegistered){
    SearchStates all_ss{s1, s2, s3};
    selector.add_searchstates(all_ss);

    traversal_info expected_i1{
        SearchStates{s1, s2},
        uniqueLoci{
            VariantLocus{5,1}, VariantLocus{7, 1}, VariantLocus{5, 2}
        }
    };

    traversal_info expected_i2{
            SearchStates{s3},
            uniqueLoci{ VariantLocus{9,1} }
    };

    new_uniqueSitePaths expected_map{
            {SitePath{5} ,  expected_i1},
            {SitePath{9}, expected_i2}
    };
    EXPECT_EQ(selector.usps, expected_map);
}

class MockRandomGenerator : public RandomGenerator{
public:
    MOCK_METHOD2(generate, uint32_t(uint32_t, uint32_t));
};

class MappingInstanceSelector_select : public ::testing::Test {
protected:
    PRG_Info prg_info;
    SearchStates ss{
        SearchState {
            SA_Interval{1, 1},
                    VariantSitePath{VariantLocus{7, 1}}
        },
        SearchState {
                SA_Interval{6, 6},
                VariantSitePath{VariantLocus{7, 2}}
        },
        SearchState {
            SA_Interval{2, 5},
                    VariantSitePath{},
                    VariantSitePath{}
        }
    };
};

TEST_F(MappingInstanceSelector_select, selectnonvariant_emptyMappingSelector){
    // Select the SearchState in invariant region of PRG
    using namespace ::testing;
    MockRandomGenerator r;
    EXPECT_CALL(r, generate(1, 2))
        .Times(Exactly(1))
        .WillOnce(Return(1));

    MappingInstanceSelector m{ss, &prg_info, &r};

    EXPECT_EQ(m.navigational_search_states.size(), 0);
}

TEST_F(MappingInstanceSelector_select, selectvariant_nonemptyMappingSelector){
    using namespace ::testing;
    MockRandomGenerator r;
    EXPECT_CALL(r, generate(1, 2))
            .Times(Exactly(1))
            .WillOnce(Return(2));

    MappingInstanceSelector m{ss, &prg_info, &r};
    EXPECT_EQ(m.navigational_search_states.size(), 2);
    uniqueLoci expected_loci{
            {VariantLocus{7, 1}},
            {VariantLocus{7, 2}}
    };
    EXPECT_EQ(m.equivalence_class_loci, expected_loci);
}
