#include <cctype>

#include "gtest/gtest.h"

#include "../../test_utils.hpp"
#include "quasimap/coverage/allele_base.hpp"
#include "quasimap/quasimap.hpp"


/*
PRG: gct5c6g6t5ag7t8c7ct
i	F	BWT	text   SA	suffix
0	0	4	3	   19	0
1	1	5	2	   10	1 3 7 4 8 2 7 2 4 0
2	2	7	4	   17	2 4 0
3	2	3	5	   1	2 4 5 2 6 3 6 4 5 1 3 7 4 8 2 7 2 4 0
4	2	5	2	   4	2 6 3 6 4 5 1 3 7 4 8 2 7 2 4 0
5	2	8	6	   15	2 7 2 4 0
6	3	0	3	   0	3 2 4 5 2 6 3 6 4 5 1 3 7 4 8 2 7 2 4 0
7	3	6	6	   6	3 6 4 5 1 3 7 4 8 2 7 2 4 0
8	3	1	4	   11	3 7 4 8 2 7 2 4 0
9	4	2	5	   18	4 0
10	4	6	1	   8	4 5 1 3 7 4 8 2 7 2 4 0
11	4	2	3	   2	4 5 2 6 3 6 4 5 1 3 7 4 8 2 7 2 4 0
12	4	7	7	   13	4 8 2 7 2 4 0
13	5	4	4	   9	5 1 3 7 4 8 2 7 2 4 0
14	5	4	8	   3	5 2 6 3 6 4 5 1 3 7 4 8 2 7 2 4 0
15	6	2	2	   5	6 3 6 4 5 1 3 7 4 8 2 7 2 4 0
16	6	3	7	   7	6 4 5 1 3 7 4 8 2 7 2 4 0
17	7	2	2	   16	7 2 4 0
18	7	3	4	   12	7 4 8 2 7 2 4 0
19	8	4	0	   14	8 2 7 2 4 0
*/

TEST(AlleleBaseCoverage, ReadCoversTwoSites_CorrectAlleleBaseCoverage) {
    auto prg_raw = "gct5c6g6t5ag7t8c7ct";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_structure(prg_info);

    uint64_t read_length = 150;

    SearchState search_state = {
            SA_Interval {3, 3},
            VariantSitePath {
                    VariantSite {5, 2},
                    VariantSite {7, 2},
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{0}, {1}, {0}},
            {{0}, {1}},
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: gct5c6g6t5ag7t8cc7ct
i	F	BWT	text	SA	suffix
0	0	4	3	    20	0
1	1	5	2	    10	1 3 7 4 8 2 2 7 2 4 0
2	2	8	4	    15	2 2 7 2 4 0
3	2	7	5	    18	2 4 0
4	2	3	2	    1	2 4 5 2 6 3 6 4 5 1 3 7 4 8 2 2 7 2 4 0
5	2	5	6	    4	2 6 3 6 4 5 1 3 7 4 8 2 2 7 2 4 0
6	2	2	3	    16	2 7 2 4 0
7	3	0	6	    0	3 2 4 5 2 6 3 6 4 5 1 3 7 4 8 2 2 7 2 4 0
8	3	6	4	    6	3 6 4 5 1 3 7 4 8 2 2 7 2 4 0
9	3	1	5	    11	3 7 4 8 2 2 7 2 4 0
10	4	2	1	    19	4 0
11	4	6	3	    8	4 5 1 3 7 4 8 2 2 7 2 4 0
12	4	2	7	    2	4 5 2 6 3 6 4 5 1 3 7 4 8 2 2 7 2 4 0
13	4	7	4	    13	4 8 2 2 7 2 4 0
14	5	4	8	    9	5 1 3 7 4 8 2 2 7 2 4 0
15	5	4	2	    3	5 2 6 3 6 4 5 1 3 7 4 8 2 2 7 2 4 0
16	6	2	2	    5	6 3 6 4 5 1 3 7 4 8 2 2 7 2 4 0
17	6	3	7	    7	6 4 5 1 3 7 4 8 2 2 7 2 4 0
18	7	2	2	    17	7 2 4 0
19	7	3	4	    12	7 4 8 2 2 7 2 4 0
20	8	4	0	    14	8 2 2 7 2 4 0
*/

TEST(AlleleBaseCoverage, ShortReadStartingOutsideSiteCoversTwoSites_FinishesBeforeSecondAlleleEnd) {
    auto prg_raw = "gct5c6g6t5ag7t8cc7ct";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_structure(prg_info);

    uint64_t read_length = 6;

    SearchState search_state = {
            SA_Interval {4, 4},
            VariantSitePath {
                    VariantSite {5, 2},
                    VariantSite {7, 2},
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
    auto prg_raw = "gct5c6g6t5ag7t8cc7ct";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_structure(prg_info);

    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {11, 11},
            VariantSitePath {
                    VariantSite {5, 3},
                    VariantSite {7, 2},
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
    auto prg_raw = "gct5c6g6t5ag7t8cc7ct";
    auto prg_info = generate_prg_info(prg_raw);

    uint64_t first_site_marker = 5;
    uint64_t second_site_marker = 7;
    auto result = inter_site_base_count(first_site_marker, second_site_marker, prg_info);
    uint64_t expected = 2;
    EXPECT_EQ(result, expected);
}

/*
PRG: ac5gg6aga5c
i	F	BWT	text	SA	suffix
0	0	2	1	    11	0
1	1	0	2	    0	1 2 5 3 3 6 1 3 1 5 2 0
2	1	6	5	    6	1 3 1 5 2 0
3	1	3	3	    8	1 5 2 0
4	2	5	3	    10	2 0
5	2	1	6	    1	2 5 3 3 6 1 3 1 5 2 0
6	3	1	1	    7	3 1 5 2 0
7	3	5	3	    3	3 3 6 1 3 1 5 2 0
8	3	3	1	    4	3 6 1 3 1 5 2 0
9	5	1	5	    9	5 2 0
10	5	2	2	    2	5 3 3 6 1 3 1 5 2 0
11	6	3	0	    5	6 1 3 1 5 2 0
*/

TEST(AlleleBaseCoverage, SaIntervalGreaterThanOne_CorrectCumulativeBaseCoverage) {
    auto prg_raw = "ac5gg6aga5c";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_structure(prg_info);

    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {7, 8},
            VariantSitePath {
                    VariantSite {5, 1},
            },
    };
    SearchStates search_states = {search_state};
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);

    auto &result = coverage.allele_base_coverage;
    SitesAlleleBaseCoverage expected = {
            {{1, 2}, {0, 0, 0}}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, ReadStartsBeforeSiteCoversFirstAllele_CorrectBaseCoverage) {
    auto prg_raw = "ac5gg6aga5c";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_structure(prg_info);

    uint64_t read_length = 150;

    SearchState search_state = {
            SA_Interval {1, 1},
            VariantSitePath {
                    VariantSite {5, 1}
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
    auto prg_raw = "ac5gg6aga5c";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_structure(prg_info);

    uint64_t read_length = 150;

    SearchState search_state = {
            SA_Interval {8, 8},
            VariantSitePath {
                    VariantSite {5, 1}
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
    auto prg_raw = "ac5gg6aga5c";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_structure(prg_info);

    uint64_t read_length = 150;

    SearchState search_state = {
            SA_Interval {6, 6},
            VariantSitePath {
                    VariantSite {5, 2}
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
    auto prg_raw = "ac5gg6aga5c";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_structure(prg_info);

    uint64_t read_length = 4;

    SearchState search_state = {
            SA_Interval {1, 1},
            VariantSitePath {
                    VariantSite {5, 2}
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
    auto prg_raw = "5gg6aga5c";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_base_structure(prg_info);
    SitesAlleleBaseCoverage expected = {
            AlleleCoverage {
                    BaseCoverage {0, 0},
                    BaseCoverage {0, 0, 0},
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenOneVariantSite_CorrectAlleleBaseCoverageStructure) {
    auto prg_raw = "ct5gg6aga5c";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_base_structure(prg_info);
    SitesAlleleBaseCoverage expected = {
            AlleleCoverage {
                    BaseCoverage {0, 0},
                    BaseCoverage {0, 0, 0},
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleBaseCoverage, GivenTwoVariantSites_CorrectAlleleBaseCoverageStructure) {
    auto prg_raw = "ct5gg6aga5ccccc7a8ttt7";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = coverage::generate::allele_base_structure(prg_info);
    SitesAlleleBaseCoverage expected = {
            AlleleCoverage {
                    BaseCoverage {0, 0},
                    BaseCoverage {0, 0, 0},
            },
            AlleleCoverage {
                    BaseCoverage {0},
                    BaseCoverage {0, 0, 0},
            },
    };
    EXPECT_EQ(result, expected);
}