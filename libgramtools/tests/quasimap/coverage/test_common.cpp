#include <cctype>
#include "gtest/gtest.h"

#include "quasimap/coverage/common.hpp"


TEST(CoverageCommon, GivenSearchStatesAndFixedSeed_CorrectRandomlySelectedSaInterval) {
    SearchStates search_states = {
            SearchState {
                    SA_Interval {1, 2}
            },
            SearchState {
                    SA_Interval {3, 4}
            }
    };
    uint32_t random_seed = 42;
    auto result = random_select_sa_interval(search_states, random_seed);
    SA_Interval expected = {1, 2};
    EXPECT_EQ(result, expected);
}


TEST(CoverageCommon, GivenTargetSaInterval_ReturnOnlySearchStatesWithTarget) {
    SearchStates search_states = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {VariantSite {11, 22}}
            },
            SearchState {
                    SA_Interval {3, 4}
            },
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {VariantSite {33, 44}}
            },
    };

    SA_Interval target_sa_interval = {1, 2};
    auto result = filter_for_sa_interval(target_sa_interval, search_states);
    SearchStates expected = {
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {VariantSite {11, 22}}
            },
            SearchState {
                    SA_Interval {1, 2},
                    VariantSitePath {VariantSite {33, 44}}
            },
    };
    EXPECT_EQ(result, expected);
}