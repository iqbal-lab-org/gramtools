#include <cctype>
#include "gtest/gtest.h"

#include "../../test_utils.hpp"
#include "quasimap/coverage/common.hpp"


using namespace gram;


/*
PRG: aa5t6cagtagcagt5ta
i	F	BTW	text	SA	suffix
0	0	1	1	    18	0
1	1	4	1	    17	1 0
2	1	0	5	    0	1 1 5 4 6 2 1 3 4 1 3 2 1 3 4 5 4 1 0
3	1	4	4	    9	1 3 2 1 3 4 5 4 1 0
4	1	2	6	    6	1 3 4 1 3 2 1 3 4 5 4 1 0
5	1	2	2	    12	1 3 4 5 4 1 0
6	1	1	1	    1	1 5 4 6 2 1 3 4 1 3 2 1 3 4 5 4 1 0
7	2	6	3	    5	2 1 3 4 1 3 2 1 3 4 5 4 1 0
8	2	3	4	    11	2 1 3 4 5 4 1 0
9	3	1	1	    10	3 2 1 3 4 5 4 1 0
10	3	1	3	    7	3 4 1 3 2 1 3 4 5 4 1 0
11	3	1	2	    13	3 4 5 4 1 0
12	4	5	1	    16	4 1 0
13	4	3	3	    8	4 1 3 2 1 3 4 5 4 1 0
14	4	3	4	    14	4 5 4 1 0
15	4	5	5	    3	4 6 2 1 3 4 1 3 2 1 3 4 5 4 1 0
16	5	4	4	    15	5 4 1 0
17	5	1	1	    2	5 4 6 2 1 3 4 1 3 2 1 3 4 5 4 1 0
18	6	4	0	    4	6 2 1 3 4 1 3 2 1 3 4 5 4 1 0
*/

TEST(CheckAlleleEncapsulated, TwoAlleleEncapsulatedMappings_True) {
    auto prg_raw = "aa5t6cagtagcagt5ta";
    auto prg_info = generate_prg_info(prg_raw);

    // read: cagt
    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {7, 8},
            VariantSitePath {
                    VariantLocus {5, 2},
            },
            SearchVariantSiteState::within_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_TRUE(result);
}


TEST(CheckAlleleEncapsulated, OneAlleleEncapsulatedMapping_True) {
    auto prg_raw = "aa5t6cagtagcagt5ta";
    auto prg_info = generate_prg_info(prg_raw);

    // read: cagt
    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {7, 7},
            VariantSitePath {
                    VariantLocus {5, 2},
            },
            SearchVariantSiteState::within_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_TRUE(result);
}


TEST(CheckAlleleEncapsulated, ReadOutsideOfSite_False) {
    auto prg_raw = "aa5t6cagtagcagt5ta";
    auto prg_info = generate_prg_info(prg_raw);

    // read: aa
    uint64_t read_length = 2;

    SearchState search_state = {
            SA_Interval {2, 2},
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_FALSE(result);
}


TEST(CheckAlleleEncapsulated, MappingExtendsOneBaseRightOustideOfSite_False) {
    auto prg_raw = "aa5t6cagtagcagt5ta";
    auto prg_info = generate_prg_info(prg_raw);

    // read: agtt
    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {5, 5},
            VariantSitePath {
                    VariantLocus {5, 2},
            },
            SearchVariantSiteState::within_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_FALSE(result);
}


TEST(CheckAlleleEncapsulated, MappingExtendsOneBaseLeftOustideOfSite_False) {
    auto prg_raw = "aa5t6cagtagcagt5ta";
    auto prg_info = generate_prg_info(prg_raw);

    // read: aca
    uint64_t read_length = 3;

    SearchState search_state = {
            SA_Interval {6, 6},
            VariantSitePath {
                    VariantLocus {5, 2},
            },
            SearchVariantSiteState::outside_variant_site
    };

    auto result = check_allele_encapsulated(search_state, read_length, prg_info);
    EXPECT_FALSE(result);
}


TEST(RrandomIntInclusive, RandomCall_MinBoundaryReturned) {
    uint64_t random_seed = 48;
    uint64_t result = random_int_inclusive(1, 10, random_seed);
    uint64_t expected = 1;
    EXPECT_EQ(result, expected);
}


TEST(RrandomIntInclusive, RandomCall_MaxBoundaryReturned) {
    uint64_t random_seed = 56;
    uint64_t result = random_int_inclusive(1, 10, random_seed);
    uint64_t expected = 10;
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
    uint64_t result = count_nonvariant_search_states(search_states);
    uint64_t expected = 1;
    EXPECT_EQ(result, expected);
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
    auto result = get_unique_path_sites(search_states);
    std::set<PathSites> expected = {
            PathSites {5, 7},
            PathSites {9, 11}
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
    auto result = get_unique_path_sites(search_states);
    std::set<PathSites> expected = {
            PathSites {9, 11}
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
    auto result = get_unique_path_sites(search_states);
    std::set<PathSites> expected = {
            PathSites {9, 11}
    };
    EXPECT_EQ(result, expected);
}


TEST(FilterForPathSites, TwoSearchStatesDifferentPaths_CorrectSingleSearchState) {
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
                            VariantLocus {9, 3},
                            VariantLocus {11, 5},
                    }
            }

    };
    PathSites target_path = {9, 11};
    auto result = filter_for_path_sites(target_path, search_states);
    SearchStates expected = {
            SearchState {
                    SA_Interval {3, 4},
                    VariantSitePath {
                            VariantLocus {9, 3},
                            VariantLocus {11, 5},
                    }
            }

    };
    EXPECT_EQ(result, expected);
}


TEST(FilterForPathSites, TwoSearchStatesDifferentPaths_CorrectEmptySearchStates) {
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
                            VariantLocus {9, 3},
                            VariantLocus {11, 5},
                    }
            }

    };
    PathSites target_path = {13, 15};
    auto result = filter_for_path_sites(target_path, search_states);
    SearchStates expected = {};
    EXPECT_EQ(result, expected);
}