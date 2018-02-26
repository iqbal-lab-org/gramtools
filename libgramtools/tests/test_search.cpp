#include <cctype>

#include "gtest/gtest.h"

#include "prg/prg.hpp"
#include "kmer_index/kmer_index.hpp"
#include "search/search.hpp"
#include "test_utils.hpp"


/*
raw PRG: gcgctggagtgctgt
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


a sa interval: 1, 2
sa_right_of_marker: 1
marker_char: 5
15, 15

sa_right_of_marker: 2
marker_char: 6
17, 17


c sa interval: 3, 7
marker_index: 7
marker_char: 5


g sa interval: 8, 11
marker_index: 11
marker_char: 6


t sa interval: 12, 14
N/A
*/


TEST(BaseNextSaInterval, GivenC_ProcessNextCharG_CorrectSaInterval) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    Marker next_char = 3;
    SA_Index next_char_first_sa_index = 8;
    SA_Interval current_sa_interval = {3, 7};

    auto result = base_next_sa_interval(next_char,
                                        next_char_first_sa_index,
                                        current_sa_interval,
                                        prg_info);
    SA_Interval expected = {8, 9};
    EXPECT_EQ(result, expected);
}


TEST(BaseNextSaInterval, GivenG_ProcessNextCharA_CorrectSaInterval) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    Marker next_char = 1;
    SA_Index next_char_first_sa_index = 1;
    SA_Interval current_sa_interval = {8, 11};

    auto result = base_next_sa_interval(next_char,
                                        next_char_first_sa_index,
                                        current_sa_interval,
                                        prg_info);
    SA_Interval expected = {1, 1};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSearch, GivenCharA_ReturnTwoCorrectSearchResults) {
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


TEST(Search, SingleCharAllele_SiteStartBoundarySingleSearchState) {
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


TEST(Search, FirstAlleleSingleChar_CorrectSkipToSiteStartBoundaryMarker) {
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


TEST(Search, CharAfterSiteEndAndAllele_FourDifferentSearchStates) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 2}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    EXPECT_EQ(markers_search_states.size(), 4);
}


TEST(Search, GivenBoundaryMarkerAndThreeAlleles_GetAlleleMarkerSaInterval) {
    auto prg_raw = "gcgct5c6g6a5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    auto boundary_marker = 5;

    auto result = get_allele_marker_sa_interval(boundary_marker, prg_info);
    SA_Interval expected = {17, 18};
    EXPECT_EQ(result, expected);
}


TEST(Search, GivenBoundaryMarkerAndTwoAlleles_GetAlleleMarkerSaInterval) {
    auto prg_raw = "aca5g6t5gcatt";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = get_allele_marker_sa_interval(5, prg_info);
    SA_Interval expected = {13, 13};
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


TEST(Search, CharAfterBoundaryEndMarker_ReturnedCorrectMarkerChars) {
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
    std::unordered_set<uint64_t> expected = {15, 17, 18};
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
        auto end_sa_index = search_state.sa_interval.first;
        auto num_chars_in_sa_interval = end_sa_index - start_sa_index + 1;
        result.push_back(num_chars_in_sa_interval);
    }
    std::vector<uint64_t> expected = {1, 1, 1};
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

    std::vector<VariantSite> result;
    for (const auto &search_state: markers_search_states)
        result.push_back(search_state.cached_variant_site);

    std::vector<VariantSite> expected = {
            {5, 1},
            {5, 2},
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
        result.push_back(search_state.cache_populated);

    std::vector<bool> expected = {true, true, true};
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


TEST(Search, ThirdAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
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
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site,
            true,
            VariantSite {5, 3}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, SecondAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
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
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site,
            true,
            VariantSite {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, FirstAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
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
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site,
            true,
            VariantSite {5, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, GivenSearchStateExitingSiteAndNextChar_CachedVariantSiteRecordedInPathHistory) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    auto pattern_char = encode_dna_base('t');

    SearchState initial_search_state = {
            SA_Interval {16, 16},
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site,
            true,
            VariantSite {5, 2}
    };
    SearchStates initial_search_states = {initial_search_state};

    auto final_search_states = search_base_backwards(pattern_char,
                                                     initial_search_states,
                                                     prg_info);

    EXPECT_EQ(final_search_states.size(), 1);
    auto search_state = final_search_states.front();
    const auto &result = search_state.variant_site_path.front();
    VariantSite expected = {5, 2};
    EXPECT_EQ(result, expected);
}


TEST(Search, InitialStateWithPopulatedVariantSitePath_CorrectVariantSitePathInResult) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);
    auto pattern_char = encode_dna_base('t');

    SearchState initial_search_state = {
            SA_Interval {16, 16},
            VariantSitePath {VariantSite {42, 53}},
            SearchVariantSiteState::outside_variant_site,
            true,
            VariantSite {5, 2}
    };
    SearchStates initial_search_states = {initial_search_state};

    auto final_search_states = search_base_backwards(pattern_char,
                                                     initial_search_states,
                                                     prg_info);

    EXPECT_EQ(final_search_states.size(), 1);
    auto search_state = final_search_states.front();
    const auto &result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantSite {5, 2},
            VariantSite {42, 53}
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


TEST(Search, GivenRead_CorrectResultSaInterval) {
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


TEST(Search, GivenReadStartingInAllele_CorrectVariantSitePath) {
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
            VariantSite {5, 3}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, GivenReadEndingInAllele_CorrectVariantSitePath) {
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
            VariantSite {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, GivenReadCrossingAllele_CorrectVariantSitePath) {
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
            VariantSite {5, 2}
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


TEST(Search, GivenReadCrossingTwoAlleles_CorrectVariantSitePath) {
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
            VariantSite {5, 1},
            VariantSite {7, 1}
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
            VariantSite {5, 1},
            VariantSite {7, 1}
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
            VariantSite {7, 2}
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
            VariantSite {7, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, KmerEndsWithinAllele_ReadCoversCorrectPath) {
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
            VariantSite {5, 3},
            VariantSite {7, 1}
    };
    EXPECT_EQ(result, expected);
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
            VariantSite {5, 3},
            VariantSite {7, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, NoMarkersToLeft_SkippingToMarkerTrue) {
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

    if (USE_SKIP_OPTIMIZATION) {
        auto result = search_state.skipping_to_marker;
        EXPECT_TRUE(result);
    }
}


TEST(Search, NoMarkersToLeft_DistanceToNextMarkerGreaterThanOne) {
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

    if (USE_SKIP_OPTIMIZATION) {
        auto result = search_state.distance_to_next_marker;
        ASSERT_TRUE(result > 1);
    }
}


/*
// TODO: remove these unit tests, they relate to SKIP_OPTIMIZATION which will
// also probably be removed
TEST(Search, MarkerToLeft_DistanceToNextMarkerIsKmerStartTextIndex) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("gtcct");
    Pattern kmer = encode_dna_bases("tcct");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    if (USE_SKIP_OPTIMIZATION) {
        auto search_state = search_states.front();
        auto result = search_state.distance_to_next_marker;
        auto expected = 2;
        EXPECT_EQ(result, expected);
    }
}


TEST(Search, MarkerToImmediateLeft_SkippingToMarkerFlagFalse) {
    auto prg_raw = "gcgct5c6g6t5agtcct";
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("agtcc");
    Pattern kmer = encode_dna_bases("gtcc");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    ASSERT_EQ(search_states.size(), 1);

    auto search_state = search_states.front();

    auto result = search_state.skipping_to_marker;
    EXPECT_FALSE(result);
}
*/


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
