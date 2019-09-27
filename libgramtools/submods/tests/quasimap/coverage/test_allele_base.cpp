#include <cctype>

#include "gtest/gtest.h"


#include "src_common/generate_prg.hpp"
#include "quasimap/coverage/allele_base.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"


using namespace gram;

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
