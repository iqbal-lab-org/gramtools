#include <cctype>
#include "gtest/gtest.h"

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
    SitePath target_path = {9, 11};
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
    SitePath target_path = {13, 15};
    auto result = filter_for_path_sites(target_path, search_states);
    SearchStates expected = {};
    EXPECT_EQ(result, expected);
}