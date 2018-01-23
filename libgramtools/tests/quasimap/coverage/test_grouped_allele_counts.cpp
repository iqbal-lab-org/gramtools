#include <cctype>
#include "gtest/gtest.h"

#include "quasimap/coverage/grouped_allele_counts.hpp"
#include "quasimap/coverage/common.hpp"
#include "../../test_utils.hpp"


TEST(GroupedAlleleCount, GivenTwoVariantSites_CorrectEmptySitesVectorSize) {
    auto prg_raw = "gct5c6g6t5ac7cc8a7";
    auto prg_info = generate_prg_info(prg_raw);
    auto grouped_allele_counts = coverage::generate::grouped_allele_counts(prg_info);

    auto result = grouped_allele_counts.size();
    uint64_t expected = 2;
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, GivenTwoSearchStates_CorrectCoverage) {
    auto prg_raw = "gct5c6g6t5ac7cc8a7";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    SearchStates search_states = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 1},
                            VariantSite {7, 1}
                    }
            },
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 2},
                            VariantSite {7, 1},
                    },
            },
    };
    coverage::record::grouped_allele_counts(coverage, search_states);
    auto result = coverage.grouped_allele_counts;
    SitesGroupedAlleleCounts expected = {
            GroupedAlleleCounts {{AlleleIds {1, 2}, 1}},
            GroupedAlleleCounts {{AlleleIds {1}, 1}},
    };
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, GivenUnorderedSearchStates_CorrectlyOrderedCoverageAlleleIds) {
    auto prg_raw = "gct5c6g6t5ac7cc8a7";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    SearchStates search_states = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 3},
                            VariantSite {7, 2}
                    }
            },
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 1},
                            VariantSite {7, 1},
                    },
            },
    };
    coverage::record::grouped_allele_counts(coverage, search_states);
    auto result = coverage.grouped_allele_counts;
    SitesGroupedAlleleCounts expected = {
            GroupedAlleleCounts {{AlleleIds {1, 3}, 1}},
            GroupedAlleleCounts {{AlleleIds {1, 2}, 1}},
    };
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, GivenSingleSearchState_CorrectCoverage) {
    auto prg_raw = "gct5c6g6t5ac7cc8a7";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    SearchStates search_states = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 3}
                    }
            }
    };
    coverage::record::grouped_allele_counts(coverage, search_states);
    auto result = coverage.grouped_allele_counts;
    SitesGroupedAlleleCounts expected = {
            GroupedAlleleCounts {{AlleleIds {3}, 1}},
            GroupedAlleleCounts {},
    };
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, MultipleSetsOfSearchStates_CorrectCoverage) {
    auto prg_raw = "gct5c6g6t5ac7cc8a7";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    SearchStates first_search_states = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 3}
                    }
            },
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 1},
                            VariantSite {7, 2}
                    },
            },
    };

    SearchStates second_search_states = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 4}
                    }
            },
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {
                            VariantSite {5, 1},
                            VariantSite {7, 2}
                    },
            },
    };

    coverage::record::grouped_allele_counts(coverage,
                                            first_search_states);
    coverage::record::grouped_allele_counts(coverage,
                                            second_search_states);

    auto result = coverage.grouped_allele_counts;
    SitesGroupedAlleleCounts expected = {
            GroupedAlleleCounts {{AlleleIds {1, 3}, 1}, {AlleleIds {1, 4}, 1}},
            GroupedAlleleCounts {{AlleleIds {2}, 2}}
    };
    EXPECT_EQ(result, expected);
}