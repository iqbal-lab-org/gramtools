/**
 * @file
 * Unit tests for vBWT backward searching.
 * Terminology:
 *  - A variant locus is where you find variant **markers**;
 *  = pairs of site & allele markers.
 *  - Search is assumed backwards; so saying we end in a site means the beginning (5' end)
 *  of the read maps there.
 *
 * Test suites:
 *  - NoVarSiteBSearch: checking regular backward searching, with no variant site markers.
 *  - MarkerSearch: checking finding and positioning variant markers in the PRG string
 *  - MarkerSAIntervals: Recovering SA Interval of variant markers.
 *  - VariantLocus_Path: checking search recovers right variant site/allele combinations.
 *  - EndInLocus: checking when search ends inside variant locus.
 *
 *  - Search: test that is not sub-classified [TODO]
 */
#include <cctype>

#include "gtest/gtest.h"

#include "prg/prg.hpp"
#include "kmer_index/build.hpp"
#include "search/search.hpp"
#include "test_utils.hpp"


using namespace gram;


/*
PRG: gcgctggagtgctgt
F -> first char of SA

i	F	BTW	text	SA
0	0	4	g	0
1	1	3	c	1 3 4 3 2 4 3 4 0
2	2	3	g	2 3 2 4 3 3 1 3 4 3 2 4 3 4 0
3	2	3	c	2 4 3 3 1 3 4 3 2 4 3 4 0
4	2	3	t	2 4 3 4 0
5	3	3	g	3 1 3 4 3 2 4 3 4 0
6	3	0	g	3 2 3 2 4 3 3 1 3 4 3 2 4 3 4 0
7	3	2	a	3 2 4 3 3 1 3 4 3 2 4 3 4 0
8	3	4	g	3 2 4 3 4 0
9	3	4	t	3 3 1 3 4 3 2 4 3 4 0
10	3	4	g	3 4 0
11	3	1	c	3 4 3 2 4 3 4 0
12	4	3	t	4 0
13	4	3	g	4 3 2 4 3 4 0
14	4	2	t	4 3 3 1 3 4 3 2 4 3 4 0
15	4	2	0	4 3 4 0
*/


TEST(Search, SingleChar_CorrectSaIntervalReturned) {
    auto prg_raw = "gcgctggagtgctgt";
    auto prg_info = generate_prg_info(prg_raw);
    auto pattern_char = encode_dna_base('g');

    SearchState initial_search_state = {
            SA_Interval {0, prg_info.fm_index.size() - 1}
    };
    SearchStates search_states = {initial_search_state};

    auto result = search_base_backwards(pattern_char,
                                        search_states,
                                        prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {5, 11},
                    VariantSitePath {}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, TwoConsecutiveChars_CorrectFinalSaIntervalReturned) {
    auto prg_raw = "gcgctggagtgctgt";
    auto prg_info = generate_prg_info(prg_raw);

    SearchState initial_search_state = {
            SA_Interval {0, prg_info.fm_index.size() - 1}
    };
    SearchStates initial_search_states = {initial_search_state};

    auto first_char = encode_dna_base('g');
    auto first_search_states = search_base_backwards(first_char,
                                                     initial_search_states,
                                                     prg_info);

    auto second_char = encode_dna_base('t');
    auto second_search_states = search_base_backwards(second_char,
                                                      first_search_states,
                                                      prg_info);

    const auto &result = second_search_states;
    SearchStates expected = {
            SearchState {
                    SA_Interval {13, 15},
                    VariantSitePath {}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, SingleCharFreqOneInText_SingleSA) {
    auto prg_raw = "gcgctggagtgctgt";
    auto prg_info = generate_prg_info(prg_raw);
    auto pattern_char = encode_dna_base('a');

    SearchState initial_search_state = {
            SA_Interval {0, prg_info.fm_index.size() - 1}
    };
    SearchStates search_states = {initial_search_state};

    auto result = search_base_backwards(pattern_char,
                                        search_states,
                                        prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {1, 1},
                    VariantSitePath {}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, TwoConsecutiveChars_SingleSaIntervalEntry) {
    auto prg_raw = "gcgctggagtgctgt";
    auto prg_info = generate_prg_info(prg_raw);

    SearchState initial_search_state = {
            SA_Interval {0, prg_info.fm_index.size() - 1}
    };
    SearchStates initial_search_states = {initial_search_state};

    auto first_char = encode_dna_base('a');
    auto first_search_states = search_base_backwards(first_char,
                                                     initial_search_states,
                                                     prg_info);

    auto second_char = encode_dna_base('g');
    auto second_search_states = search_base_backwards(second_char,
                                                      first_search_states,
                                                      prg_info);

    const auto &result = second_search_states.front().sa_interval;
    SA_Interval expected{5, 5};
    EXPECT_EQ(result, expected);
}


TEST(Search, TwoConsecutiveCharsNoValidSaInterval_NoSearchStatesReturned) {
    auto prg_raw = "gcgctggagtgctgt";
    auto prg_info = generate_prg_info(prg_raw);

    SearchState initial_search_state = {
            SA_Interval {0, prg_info.fm_index.size() - 1}
    };
    SearchStates initial_search_states = {initial_search_state};

    auto first_char = encode_dna_base('a');
    auto first_search_states = search_base_backwards(first_char,
                                                     initial_search_states,
                                                     prg_info);

    auto second_char = encode_dna_base('c');
    const auto &result = search_base_backwards(second_char,
                                               first_search_states,
                                               prg_info);

    SearchStates expected = {};
    EXPECT_EQ(result, expected);
}


/*
PRG: gcgct5c6g6a5agtcct

i   F   BTW text  SA   suffix
0   0   4   3     18     0
1   1   5   2     12     1 3 4 2 2 4 0
2   1   6   3     10     1 5 1 3 4 2 2 4 0
3   2   4   2     15     2 2 4 0
4   2   3   4     1      2 3 2 4 5 2 6 3 6 1 5 1 3 4 2 2 4 0
5   2   2   5     16     2 4 0
6   2   3   2     3      2 4 5 2 6 3 6 1 5 1 3 4 2 2 4 0
7   2   5   6     6      2 6 3 6 1 5 1 3 4 2 2 4 0
8   3   0   3     0      3 2 3 2 4 5 2 6 3 6 1 5 1 3 4 2 2 4 0
9   3   2   6     2      3 2 4 5 2 6 3 6 1 5 1 3 4 2 2 4 0
10  3   1   1     13     3 4 2 2 4 0
11  3   6   5     8      3 6 1 5 1 3 4 2 2 4 0
12  4   2   1     17     4 0
13  4   3   3     14     4 2 2 4 0
14  4   2   4     4      4 5 2 6 3 6 1 5 1 3 4 2 2 4 0
15  5   1   2     11     5 1 3 4 2 2 4 0
16  5   4   2     5      5 2 6 3 6 1 5 1 3 4 2 2 4 0
17  6   3   4     9      6 1 5 1 3 4 2 2 4 0
18  6   2   0     7      6 3 6 1 5 1 3 4 2 2 4 0

*/


TEST(NoVarSiteBSearch, GivenC_ProcessNextCharG_CorrectSaInterval) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    Marker next_char = 3;
    SA_Index next_char_first_sa_index = 8;
    SA_Interval current_sa_interval = {3, 7}; // all C

    auto result = base_next_sa_interval(next_char,
                                        next_char_first_sa_index,
                                        current_sa_interval,
                                        prg_info);
    SA_Interval expected = {8, 9};
    EXPECT_EQ(result, expected);
}


TEST(NoVarSiteBSearch, GivenG_ProcessNextCharA_CorrectSaInterval) {
    // Looking for 'ag' here
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    Marker next_char = 1;
    SA_Index next_char_first_sa_index = 1;
    SA_Interval current_sa_interval = {8, 11}; // all G

    auto result = base_next_sa_interval(next_char,
                                        next_char_first_sa_index,
                                        current_sa_interval,
                                        prg_info);
    SA_Interval expected = {1, 1};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSearch, GivenCharA_FindLeftMarkers_AndSeedSearchStates){
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 2}
    };

    auto result = left_markers_search(initial_search_state,
                                      prg_info);
    MarkersSearchResults expected = {
            {1, 5},
            {2, 6},
    };
    EXPECT_EQ(result, expected);

    // Expect three: one for exiting the site; two for entering.
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    EXPECT_EQ(markers_search_states.size(), 3);
}

TEST(MarkerSearch, TestSiteMarker_Entry_or_Exit){
    auto prg_raw = "gcgct5C6g6a5Agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    const Marker marker_char = 5;

    // TEST 1: char a at site exit point
    SA_Index sa_right_of_marker = 1;

    auto siteInfo = gram::site_boundary_marker_info(marker_char, sa_right_of_marker, prg_info);
    EXPECT_EQ(false, siteInfo.is_start_boundary);
    EXPECT_EQ(15, siteInfo.sa_interval.first);

    // TEST 2: char c at site entry point
    sa_right_of_marker = 7;
    siteInfo = gram::site_boundary_marker_info(marker_char, sa_right_of_marker, prg_info);
    EXPECT_EQ(true, siteInfo.is_start_boundary);
    EXPECT_EQ(16, siteInfo.sa_interval.first);
}


TEST(MarkerSearch, GivenCharG_ReturnOneCorrectSearchResults) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    // first char: g
    SearchState initial_search_state = {
            SA_Interval {8, 11}
    };

    auto result = left_markers_search(initial_search_state,
                                      prg_info);
    MarkersSearchResults expected = {
            {11, 6},
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, SingleCharAllele_CorrectSkipToSiteStartBoundaryMarker) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    // first char: g
    SearchState initial_search_state = {
            SA_Interval {8, 11}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    const auto &first_markers_search_state = markers_search_states.front();

    const auto &result = first_markers_search_state.sa_interval;
    SA_Interval expected = {16, 16};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSearch, GivenCharG_NoMarkersToLeft){
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    // first char: g
    SearchState initial_search_state = {
            SA_Interval {8, 11}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    const auto &result = markers_search_states.size();
    auto expected = 1;
    EXPECT_EQ(result, expected);
}


TEST(MarkerSearch, GivenCharC_GoToVarSiteStart) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    // first char: c
    SearchState initial_search_state = {
            SA_Interval {3, 7}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    const auto &first_markers_search_state = markers_search_states.front();

    EXPECT_EQ(markers_search_states.size(), 1);
    auto &result = first_markers_search_state.sa_interval;
    SA_Interval expected = {16, 16};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSAIntervals, BoundaryMarkerAndThreeAlleles_GetAlleleMarkerSaInterval) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    Marker boundary_marker = 5;

    auto result = get_allele_marker_sa_interval(boundary_marker, prg_info);
    SA_Interval expected = {17, 18};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSAIntervals, BoundaryMarkerAndTwoAlleles_GetAlleleMarkerSaInterval) {
    auto prg_raw = "aca5g6t5gcatt";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = get_allele_marker_sa_interval(5, prg_info);
    SA_Interval expected = {13, 13};
    EXPECT_EQ(result, expected);
}


/*
PRG: 7g8c7g9t10a9
i	F	BTW	text	SA	suffix
0	0	9	7	    11	0
1	1	10	3	    9	1 9 0
2	2	8	8	    3	2 7 3 9 4 10 1 9 0
3	3	7	2	    1	3 8 2 7 3 9 4 10 1 9 0
4	3	7	7	    7	3 9 4 10 1 9 0
5	4	9	3	    7	4 10 1 9 0
6	7	0	9	    0	7 3 8 2 7 3 9 4 10 1 9 0
7	7	2	4	    4	7 3 9 4 10 1 9 0
8	8	3	10	    2	8 2 7 3 9 4 10 1 9 0
9	9	1	1	    10	9 0
10	9	3	9	    8	9 4 10 1 9 0
11	10	4	0	    8	10 1 9 0
*/
TEST(MarkerSAIntervals, GivenPrgWithNonContinuousAlphabet_CorrectAlleleMarkerEndBoundary) {
    auto prg_raw = "7g8c7g9t10a9";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = get_allele_marker_sa_interval(7, prg_info);
    SA_Interval expected = {8, 8};
    EXPECT_EQ(result, expected);
}


/*
PRG: gcgct5c6g6t5agtcct
i	F	BWT	text	SA	suffix
0	0	4	 3	    18	  0
1	1	5	 2	    12	  1 3 4 2 2 4 0
2	2	4	 3	    15	  2 2 4 0
3	2	3	 2	    1	  2 3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
4	2	2	 4	    16	  2 4 0
5	2	3	 5	    3	  2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
6	2	5	 2	    6	  2 6 3 6 4 5 1 3 4 2 2 4 0
7	3	0	 6	    0	  3 2 3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
8	3	2	 3	    2	  3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
9	3	1	 6	    13	  3 4 2 2 4 0
10	3	6	 4	    8	  3 6 4 5 1 3 4 2 2 4 0
11	4	2	 5	    17	  4 0
12	4	3	 1	    14	  4 2 2 4 0
13	4	6	 3	    10	  4 5 1 3 4 2 2 4 0
14	4	2	 4	    4	  4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
15	5	4	 2	    11	  5 1 3 4 2 2 4 0
16	5	4	 2	    5	  5 2 6 3 6 4 5 1 3 4 2 2 4 0
17	6	2	 4	    7	  6 3 6 4 5 1 3 4 2 2 4 0
18	6	3	 0	    9	  6 4 5 1 3 4 2 2 4 0
 */


TEST(MarkerSearch, AtSiteEnd_GetAllMarkerChars) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 1}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);

    std::unordered_set<uint64_t> result;
    for (const auto &search_state: markers_search_states) {
        auto sa_index = search_state.sa_interval.first;
        auto text_index = prg_info.fm_index[sa_index];
        auto marker_char = prg_info.fm_index.text[text_index];
        result.insert(marker_char);
    }
    std::unordered_set<uint64_t> expected = {6, 6, 5};
    EXPECT_EQ(result, expected);
}


TEST(Search, CharAfterBoundaryEndMarker_ReturnedCorrectSaIndexes) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 1}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);

    std::unordered_set<uint64_t> result;
    for (const auto &search_state: markers_search_states) {
        auto start_sa_index = search_state.sa_interval.first;
        result.insert(start_sa_index);
    }
    std::unordered_set<uint64_t> expected = {15, 17};
    EXPECT_EQ(result, expected);
}


TEST(Search, CharAfterBoundaryEndMarker_ReturnedSingleCharSaIntervals) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 1}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);

    std::vector<uint64_t> result;
    for (const auto &search_state: markers_search_states) {
        auto start_sa_index = search_state.sa_interval.first;
        auto end_sa_index = search_state.sa_interval.second;
        auto num_chars_in_sa_interval = end_sa_index - start_sa_index + 1;
        result.push_back(num_chars_in_sa_interval);
    }
    std::vector<uint64_t> expected = {2, 1};
    EXPECT_EQ(result, expected);
}


TEST(Search, CharAfterBoundaryEndMarker_ReturnedSearchStatesHaveCorrectLastVariantSiteAttributes) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 1}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);

    std::vector<VariantLocus> result;
    for (const auto &search_state: markers_search_states)
        result.push_back(search_state.variant_site_path.front());

    // We expect the following: one searchstate has two alleles, and thus the allele part is unspecified still
    // The other is a singleton SearchState corresponding to the end of the site, which is alleles #3.
    std::vector<VariantLocus> expected = {
            {5, ALLELE_UNKNOWN},
            {5, 3},
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, CharAfterBoundaryEndMarker_ReturnedSearchStatesHaveCorrectVariantSiteRecordedAttributes) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 1}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    std::vector<bool> result;
    for (const auto &search_state: markers_search_states)
        result.push_back(search_state.variant_site_path.size()==1);

    std::vector<bool> expected = {true, true};
    EXPECT_EQ(result, expected);
}

TEST(Search, GivenAlleleMarkerSaIndex_ReturnAlleleId) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    uint64_t allele_marker_sa_index = 18;
    auto result = get_allele_id(allele_marker_sa_index,
                                prg_info);
    auto expected = 2;
    EXPECT_EQ(result, expected);
}


TEST(ExitASite, ThirdAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    // first char: t
    SearchState initial_search_state = {
            SA_Interval {11, 14}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    EXPECT_EQ(markers_search_states.size(), 1);
    auto result = markers_search_states.front();
    SearchState expected = {
            SA_Interval {16, 16},
            VariantSitePath {VariantLocus{5,3}},
            SearchVariantSiteState::outside_variant_site,
    };
    EXPECT_EQ(result, expected);
}


TEST(ExitASite, SecondAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    // first char: g
    SearchState initial_search_state = {
            SA_Interval {7, 10}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    EXPECT_EQ(markers_search_states.size(), 1);
    auto result = markers_search_states.front();
    SearchState expected = {
            SA_Interval {16, 16},
            VariantSitePath {VariantLocus{5, 2}},
            SearchVariantSiteState::outside_variant_site,
    };
    EXPECT_EQ(result, expected);
}


TEST(ExitASite, FirstAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    // first char: c
    SearchState initial_search_state = {
            SA_Interval {2, 6}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    EXPECT_EQ(markers_search_states.size(), 1);
    auto result = markers_search_states.front();
    SearchState expected = {
            SA_Interval {16, 16},
            VariantSitePath {VariantLocus{5,1}},
            SearchVariantSiteState::outside_variant_site,
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: gcgct5c6g6t5agtcct
i	F	BWT	text	SA	suffix
0	0	4	 3	    18	  0
1	1	5	 2	    12	  1 3 4 2 2 4 0
2	2	4	 3	    15	  2 2 4 0
3	2	3	 2	    1	  2 3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
4	2	2	 4	    16	  2 4 0
5	2	3	 5	    3	  2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
6	2	5	 2	    6	  2 6 3 6 4 5 1 3 4 2 2 4 0
7	3	0	 6	    0	  3 2 3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
8	3	2	 3	    2	  3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
9	3	1	 6	    13	  3 4 2 2 4 0
10	3	6	 4	    8	  3 6 4 5 1 3 4 2 2 4 0
11	4	2	 5	    17	  4 0
12	4	3	 1	    14	  4 2 2 4 0
13	4	6	 3	    10	  4 5 1 3 4 2 2 4 0
14	4	2	 4	    4	  4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
15	5	4	 2	    11	  5 1 3 4 2 2 4 0
16	5	4	 2	    5	  5 2 6 3 6 4 5 1 3 4 2 2 4 0
17	6	2	 4	    7	  6 3 6 4 5 1 3 4 2 2 4 0
18	6	3	 0	    9	  6 4 5 1 3 4 2 2 4 0
 */
TEST(Search, InitialStateWithPopulatedVariantSitePath_CorrectVariantSitePathInResult) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    auto pattern_char = encode_dna_base('t');

    SearchState initial_search_state = {
            SA_Interval {10,10}, //Starting at char 'c' at index 6 in prg
            VariantSitePath {},
            SearchVariantSiteState::unknown,
    };
    SearchStates initial_search_states = {initial_search_state};

    auto final_search_states=process_read_char_search_states(pattern_char,initial_search_states,prg_info);

    EXPECT_EQ(final_search_states.size(), 1);
    auto search_state = final_search_states.front();
    const auto &result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 2},
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, KmerAbsentFromKmerIndex_NoSearchStatesReturned) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("tagtaa");
    Pattern kmer = encode_dna_bases("gtaa");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 0);
}


TEST(SA_Interval, GivenRead_CorrectResultSaInterval) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("tagtcc");
    Pattern kmer = encode_dna_bases("gtcc");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    ASSERT_EQ(search_states.size(), 1);

    auto search_state = search_states.front();
    auto result = search_state.sa_interval;
    SA_Interval expected = {13, 13};
    EXPECT_EQ(result, expected);
}


TEST(VariantLocus_Path, GivenSearchEndingInAllele_CorrectVariantSitePath) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("tagtcc");
    Pattern kmer = encode_dna_bases("gtcc");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 3}
    };
    EXPECT_EQ(result, expected);
}


TEST(VariantLocus_Path, GivenSearchStartingInAllele_CorrectVariantSitePath) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("cgctg");
    Pattern kmer = encode_dna_bases("gctg");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(VariantLocus_Path, GivenSearchCrossingAllele_CorrectVariantSitePath) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("ctgag");
    Pattern kmer = encode_dna_bases("tgag");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);
}


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


TEST(VariantLocus_Path, GivenReadCrossingTwoAlleles_CorrectVariantSitePath) {
    auto prg_raw = "gct5c6g6t5ag7t8c7ct";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("tct");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("cagtct");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 1},
            VariantLocus {7, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, KmerWithinAlleleNotCrossingMarker_ReadCoversCorrectPath) {
    auto prg_raw = "gct5c6g6t5ag7tct8c7ct";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("tct");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("cagtct");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 1},
            VariantLocus {7, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, KmerImmediatelyAfterVariantSite_ReadCoversCorrectPath) {
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("cta");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("gccta");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {7, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, KmerCrossesVariantSite_ReadCoversCorrectPath) {
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("gccta");
    Patterns kmers = {kmer};
    auto kmer_size = 5;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("agccta");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {7, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(EndInLocus, SearchStarts_andEnds_withinLoci) {
    auto prg_raw = "gct5c6g6T5AG7T8c7cta";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("tagt");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 3},
            VariantLocus {7, 1}
    };
    EXPECT_EQ(result, expected);
}

TEST(EndInLocus, SearchEnds_AtOneAlleleMarker) {
    auto prg_raw = "gct5c6G6t5AG7T8c7cta";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("gagt");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 2},
            VariantLocus {7, 1}
    };
    EXPECT_EQ(result, expected);
}

/*
 * A case where we end the read mapping inside several alleles of the same site.
 * We test expected behaviour along the way from kmer indexing to read mapping alleles concurrently
 * to allele ID specification post mapping.
 */
TEST(EndInLocus, SearchEnds_AtConcurrentAlleles) {
    auto prg_raw = "gct5gC6aC6C6t5Cg";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("c");
    Patterns kmers = {kmer};
    auto kmer_size = 1;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    // KMER INDEXING
    // We expect five occurrences of 'C' at this stage, in a single SA interval
    auto search_states = kmer_index.at(kmer);
    EXPECT_EQ(search_states.size(), 1);
    SA_Interval sa = search_states.front().sa_interval;
    EXPECT_EQ(sa.second - sa.first + 1, 5);

    // Next up, look for a C
    Base pattern_char = 2;
    search_states = process_read_char_search_states(pattern_char,
                                                        search_states,
                                                        prg_info);

    // CONCURRENT ALLELE QUERYING
    // We expect three occurrences of 'CC' at this stage, in a single SA interval - because
    // the allele markers sort together in the SA. The allele IDs should be unspecified.
    EXPECT_EQ(search_states.size(), 1);
    EXPECT_EQ(search_states.front().variant_site_path.front().second, ALLELE_UNKNOWN);

    // ALLELE ID SPECIFICATION
    // This function gets called when we have finished mapping our read and we have unknown allele ids left.
    gram::set_allele_ids(search_states, prg_info);
    EXPECT_EQ(search_states.size(), 3);

    for (const auto search_state: search_states){
        SA_Interval sa = search_states.front().sa_interval;
        EXPECT_EQ(sa.second - sa.first + 1, 1);
    }
}


TEST(Search, KmerCrossesMultipleVariantSites_ReadCoversCorrectPath) {
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("tagt");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("cttagt");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantLocus {5, 3},
            VariantLocus {7, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, OneMappingEncapsulatedByAllele_StateIsWithinVariantSite) {
    auto prg_raw = "t5c6gcttagt5aa";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("tagt");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("cttagt");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_state;
    SearchVariantSiteState expected = SearchVariantSiteState::within_variant_site;
    EXPECT_EQ(result, expected);
}


TEST(Search, TwoMappingsEncapsulatedByAllele_StateIsWithinVariantSite) {
    auto prg_raw = "t5c6gcttagtacgcttagt5aa";
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("tagt");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("cttagt");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.variant_site_state;
    SearchVariantSiteState expected = SearchVariantSiteState::within_variant_site;
    EXPECT_EQ(result, expected);
}


/*
PRG: ac5t6cagtagtc5ta
i	F	BWT	text	SA	suffix
0	0	1	1	    16	0
1	1	4	2	    15	1 0
2	1	0	5	    0	1 2 5 4 6 2 1 3 4 1 3 4 2 5 4 1 0
3	1	2	4	    6	1 3 4 1 3 4 2 5 4 1 0
4	1	4	6	    9	1 3 4 2 5 4 1 0
5	2	6	2	    5	2 1 3 4 1 3 4 2 5 4 1 0
6	2	4	1	    12	2 5 4 1 0
7	2	1	3	    1	2 5 4 6 2 1 3 4 1 3 4 2 5 4 1 0
8	3	1	4	    7	3 4 1 3 4 2 5 4 1 0
9	3	1	1	    10	3 4 2 5 4 1 0
10	4	5	3	    14	4 1 0
11	4	3	4	    8	4 1 3 4 2 5 4 1 0
12	4	3	2	    11	4 2 5 4 1 0
13	4	5	5	    3	4 6 2 1 3 4 1 3 4 2 5 4 1 0
14	5	2	4	    13	5 4 1 0
15	5	2	1	    2	5 4 6 2 1 3 4 1 3 4 2 5 4 1 0
16	6	4	0	    4	6 2 1 3 4 1 3 4 2 5 4 1 0
*/
TEST(HandleAlleleEncapsulatedStates, AlleleEncapsulatedStateMissingPath_CorrectPathSet) {
    auto prg_raw = "ac5t6cagtagtc5ta";
    auto prg_info = generate_prg_info(prg_raw);
    SearchStates search_states = {
            SearchState {
                    SA_Interval {8, 8}
            }
    };
    auto result = handle_allele_encapsulated_states(search_states, prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {8, 8},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(HandleAlleleEncapsulatedStates, AlleleEncapsulatedState_NoChange) {
    auto prg_raw = "ac5t6cagtagtc5ta";
    auto prg_info = generate_prg_info(prg_raw);
    SearchStates search_states = {
            SearchState {
                    SA_Interval {8, 8},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    SearchVariantSiteState::within_variant_site
            }
    };
    auto result = handle_allele_encapsulated_states(search_states, prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {8, 8},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(HandleAlleleEncapsulatedStates, SaIntervalGreaterThanOneAlleleEncapsulated_CorrectPathSet) {
    auto prg_raw = "ac5t6cagtagtc5ta";
    auto prg_info = generate_prg_info(prg_raw);
    SearchStates search_states = {
            SearchState {
                    SA_Interval {3, 4}
            }
    };
    auto result = handle_allele_encapsulated_states(search_states, prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {3, 4},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: gcgct5c6g6t5agtcct
i	F	BWT	text	SA	suffix
0	0	4	 3	    18	  0
1	1	5	 2	    12	  1 3 4 2 2 4 0
2	2	4	 3	    15	  2 2 4 0
3	2	3	 2	    1	  2 3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
4	2	2	 4	    16	  2 4 0
5	2	3	 5	    3	  2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
6	2	5	 2	    6	  2 6 3 6 4 5 1 3 4 2 2 4 0
7	3	0	 6	    0	  3 2 3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
8	3	2	 3	    2	  3 2 4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
9	3	1	 6	    13	  3 4 2 2 4 0
10	3	6	 4	    8	  3 6 4 5 1 3 4 2 2 4 0
11	4	2	 5	    17	  4 0
12	4	3	 1	    14	  4 2 2 4 0
13	4	6	 3	    10	  4 5 1 3 4 2 2 4 0
14	4	2	 4	    4	  4 5 2 6 3 6 4 5 1 3 4 2 2 4 0
15	5	4	 2	    11	  5 1 3 4 2 2 4 0
16	5	4	 2	    5	  5 2 6 3 6 4 5 1 3 4 2 2 4 0
17	6	2	 4	    7	  6 3 6 4 5 1 3 4 2 2 4 0
18	6	3	 0	    9	  6 4 5 1 3 4 2 2 4 0
 */


TEST(HandleAlleleEncapsulatedStates, OutsideSite_NoPathSet) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    SearchStates search_states = {
            SearchState {
                    SA_Interval {7, 7}
            }
    };
    auto result = handle_allele_encapsulated_states(search_states, prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {7, 7},
                    VariantSitePath {},
                    SearchVariantSiteState::outside_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}

/*
PRG: cagtaa5t6cagtaggc5ta
i	F	BTW	text	SA	suffix
0	0	1	2	    20	0
1	1	4	1	    19	1 0
2	1	4	3	    4	1 1 5 4 6 2 1 3 4 1 3 3 2 5 4 1 0
3	1	4	4	    13	1 3 3 2 5 4 1 0
4	1	2	1	    1	1 3 4 1 1 5 4 6 2 1 3 4 1 3 3 2 5 4 1 0
5	1	2	1	    10	1 3 4 1 3 3 2 5 4 1 0
6	1	1	5	    5	1 5 4 6 2 1 3 4 1 3 3 2 5 4 1 0
7	2	0	4	    0	2 1 3 4 1 1 5 4 6 2 1 3 4 1 3 3 2 5 4 1 0
8	2	6	6	    9	2 1 3 4 1 3 3 2 5 4 1 0
9	2	3	2	    16	2 5 4 1 0
10	3	3	1	    15	3 2 5 4 1 0
11	3	1	3	    14	3 3 2 5 4 1 0
12	3	1	4	    2	3 4 1 1 5 4 6 2 1 3 4 1 3 3 2 5 4 1 0
13	3	1	1	    11	3 4 1 3 3 2 5 4 1 0
14	4	5	3	    18	4 1 0
15	4	3	3	    3	4 1 1 5 4 6 2 1 3 4 1 3 3 2 5 4 1 0
16	4	3	2	    12	4 1 3 3 2 5 4 1 0
17	4	5	5	    7	4 6 2 1 3 4 1 3 3 2 5 4 1 0
18	5	2	4	    17	5 4 1 0
19	5	1	1	    6	5 4 6 2 1 3 4 1 3 3 2 5 4 1 0
20	6	4	0	    8	6 2 1 3 4 1 3 3 2 5 4 1 0
*/

TEST(HandleAlleleEncapsulatedState, ReadAlleleEncapsulatedAndOutsideSite_SplitIntoTwoSearchStates) {
    auto prg_raw = "cagtaa5t6cagtaggc5ta";
    auto prg_info = generate_prg_info(prg_raw);

    SearchState search_state = {
            SA_Interval {7, 8}

    };
    auto result = handle_allele_encapsulated_state(search_state, prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {7, 7},
                    VariantSitePath {},
                    SearchVariantSiteState::outside_variant_site
            },
            SearchState {
                    SA_Interval {8, 8},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: tcagtt5tcagtcag6atcagtttcag5ta7atcagt8gtg7g
i	F	BWT	text	SA	suffix
0	0	3	4	    43	  0
1	1	2	2	    9	  1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
2	1	2	1	    19	  1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
3	1	2	3	    2	  1 3 4 4 5 4 2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
4	1	2	4	    34	  1 3 4 8 3 4 3 7 3 0
5	1	2	4	    25	  1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
6	1	2	5	    13	  1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
7	1	6	4	    16	  1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
8	1	7	2	    31	  1 4 2 1 3 4 8 3 4 3 7 3 0
9	1	4	1	    29	  1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
10	2	4	3	    8	  2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
11	2	4	4	    18	  2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
12	2	4	2	    1	  2 1 3 4 4 5 4 2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
13	2	4	1	    33	  2 1 3 4 8 3 4 3 7 3 0
14	2	4	3	    24	  2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
15	2	4	6	    12	  2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
16	3	7	1	    42	  3 0
17	3	1	4	    10	  3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
18	3	8	2	    38	  3 4 3 7 3 0
19	3	1	1	    20	  3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
20	3	1	3	    3	  3 4 4 5 4 2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
21	3	1	4	    35	  3 4 8 3 4 3 7 3 0
22	3	1	4	    26	  3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
23	3	1	4	    14	  3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
24	3	4	2	    40	  3 7 3 0
25	4	5	1	    28	  4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
26	4	5	3	    7	  4 2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
27	4	1	5	    17	  4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
28	4	0	4	    0	  4 2 1 3 4 4 5 4 2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
29	4	1	1	    32	  4 2 1 3 4 8 3 4 3 7 3 0
30	4	4	7	    23	  4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
31	4	3	1	    11	  4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
32	4	3	4	    39	  4 3 7 3 0
33	4	4	2	    22	  4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
34	4	3	1	    21	  4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
35	4	3	3	    4	  4 4 5 4 2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
36	4	4	4	    5	  4 5 4 2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
37	4	3	8	    36	  4 8 3 4 3 7 3 0
38	5	3	3	    27	  5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
39	5	4	4	    6	  5 4 2 1 3 4 2 1 3 6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
40	6	3	3	    15	  6 1 4 2 1 3 4 4 4 2 1 3 5 4 1 7 1 4 2 1 3 4 8 3 4 3 7 3 0
41	7	1	7	    30	  7 1 4 2 1 3 4 8 3 4 3 7 3 0
42	7	3	3	    41	  7 3 0
43	8	4	0	    37	  8 3 4 3 7 3 0
*/
TEST(HandleAlleleEncapsulatedState, MappingMultipleAlleleEncapsulation_CorrectSearchStates) {
    auto prg_raw = "tcagtt5tcagtcag6atcagtttcag5ta7atcagt8gtg7g";
    auto prg_info = generate_prg_info(prg_raw);

    SearchState search_state = {
            SA_Interval {10, 15}

    };
    auto result = handle_allele_encapsulated_state(search_state, prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {10, 10},
                    VariantSitePath {
                            VariantLocus {5, 1}
                    },
                    SearchVariantSiteState::within_variant_site
            },
            SearchState {
                    SA_Interval {11, 11},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    SearchVariantSiteState::within_variant_site
            },
            SearchState {
                    SA_Interval {12, 12},
                    VariantSitePath {},
                    SearchVariantSiteState::outside_variant_site
            },
            SearchState {
                    SA_Interval {13, 13},
                    VariantSitePath {
                            VariantLocus {7, 1}
                    },
                    SearchVariantSiteState::within_variant_site
            },
            SearchState {
                    SA_Interval {14, 14},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    SearchVariantSiteState::within_variant_site
            },
            SearchState {
                    SA_Interval {15, 15},
                    VariantSitePath {
                            VariantLocus {5, 1}
                    },
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, ReadLeadsToPrgEdge_NoSearchStatesFound) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("agcgc");
    Pattern kmer = encode_dna_bases("gcgc");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    ASSERT_TRUE(search_states.empty());
}
