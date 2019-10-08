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
 *  - VarSiteBSearch: backward searching with var site markers.
 *
 *  - MarkerSearch: checking finding and positioning variant markers in the PRG string
 *  - MarkerSAIntervals: Recovering SA Interval of variant markers.
 *
 *  - VariantLocus_Path: checking search recovers right variant site/allele combinations.
 *  - EndInLocus: checking when search ends inside variant locus.
 *  - StartEndInLocus: search starts and end inside VariantLocus.
 *
 *  - Search: test that is not sub-classified [TODO]
 */
#include <cctype>

#include "gtest/gtest.h"

#include "src_common/generate_prg.hpp"
#include "prg/prg.hpp"
#include "kmer_index/build.hpp"
#include "search/search.hpp"


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
    auto prg_raw = encode_prg("gcgctggagtgctgt");
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
    auto prg_raw = encode_prg("gcgctggagtgctgt");
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
    auto prg_raw = encode_prg("gcgctggagtgctgt");
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
    auto prg_raw = encode_prg("gcgctggagtgctgt");
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
    auto prg_raw = encode_prg("gcgctggagtgctgt");
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
PRG: GCGCT5C6G6A6AGTCCT
i	BWT	SA	text_suffix
0	G	18
1	6	12	A G T C C T
2	6	10	A 6 A G T C C T
3	G	15	C C T
4	T	1	C G C T 5 C 6 G 6 A 6 A G T C C T
5	C	16	C T
6	T	3	C T 5 C 6 G 6 A 6 A G T C C T
7	5	6	C 6 G 6 A 6 A G T C C T
8	0	0	G C G C T 5 C 6 G 6 A 6 A G T C C T
9	C	2	G C T 5 C 6 G 6 A 6 A G T C C T
10	A	13	G T C C T
11	6	8	G 6 A 6 A G T C C T
12	C	17	T
13	T	14	T C C T
14	C	4	T 5 C 6 G 6 A 6 A G T C C T
15	G	5	5 C 6 G 6 A 6 A G T C C T
16	A	11	6 A G T C C T
17	T	9	6 A 6 A G T C C T
18	C	7	6 G 6 A 6 A G T C C T
*/


TEST(NoVarSiteBSearch, GivenC_ProcessNextCharG_CorrectSaInterval) {
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
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
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
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
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
    auto prg_info = generate_prg_info(prg_raw);
    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 2}
    };

    auto result = left_markers_search(initial_search_state,
                                      prg_info);
    MarkersSearchResults expected = {
            {6, 0},
            {5, 3},
    };
    EXPECT_EQ(result, expected);

    // Expect three: one for exiting the site; two for entering.
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    EXPECT_EQ(markers_search_states.size(), 2);
}

TEST(MarkerSearch, TestSiteMarker_Entry_or_Exit){
    auto prg_raw = encode_prg("gcgct5C6g6a6Agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    const Marker marker_char = 6;

    // TEST 1: char a at site exit point
    SA_Index sa_right_of_marker = 1;
    auto text_index = prg_info.fm_index[sa_right_of_marker];

    auto isSiteEnd = gram::marker_is_site_end(marker_char, sa_right_of_marker, prg_info);
    EXPECT_EQ(true, isSiteEnd);

    // TEST 2: char c at site entry point
    sa_right_of_marker = 7;
    isSiteEnd = gram::marker_is_site_end(marker_char, sa_right_of_marker, prg_info);
    EXPECT_EQ(false, isSiteEnd);
}


TEST(MarkerSearch, GivenCharG_ReturnOneCorrectSearchResults) {
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
    auto prg_info = generate_prg_info(prg_raw);
    // first char: g
    SearchState initial_search_state = {
            SA_Interval {8, 11}
    };

    auto result = left_markers_search(initial_search_state,
                                      prg_info);
    MarkersSearchResults expected = {
            {5, 2},
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, SingleCharAllele_CorrectSkipToSiteStartBoundaryMarker) {
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
    auto prg_info = generate_prg_info(prg_raw);
    // first char: g
    SearchState initial_search_state = {
            SA_Interval {8, 11}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    const auto &first_markers_search_state = markers_search_states.front();

    const auto &result = first_markers_search_state.sa_interval;
    SA_Interval expected = {15, 15};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSearch, GivenCharG_NoMarkersToLeft){
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
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
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
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
    SA_Interval expected = {15, 15};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSAIntervals, BoundaryMarkerAndThreeAlleles_GetAlleleMarkerSaInterval) {
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
    auto prg_info = generate_prg_info(prg_raw);
    Marker allele_marker = 6;

    auto result = get_allele_marker_sa_interval(allele_marker, prg_info);
    SA_Interval expected = {16, 18};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSAIntervals, BoundaryMarkerAndTwoAlleles_GetAlleleMarkerSaInterval) {
    auto prg_raw = encode_prg("aca5g6t6catt");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = get_allele_marker_sa_interval(6, prg_info);
    SA_Interval expected = {11, 12};
    EXPECT_EQ(result, expected);
}


/*
PRG: 7G8C8G9T10A10
i	BWT	SA	text_suffix
0	10	11	1
1	10	9	0 A 1
2	8	3	C 8 G 9 T 1 0 A 1
3	7	1	G 8 C 8 G 9 T 1 0 A 1
4	8	5	G 9 T 1 0 A 1
5	9	7	T 1 0 A 1
6	0	0	7 G 8 C 8 G 9 T 1 0 A 1
7	G	2	8 C 8 G 9 T 1 0 A 1
8	C	4	8 G 9 T 1 0 A 1
9	G	6	9 T 1 0 A 1
10	A	10	A 1
11	T	8	1 0 A 1
*/
TEST(MarkerSAIntervals, GivenPrgWithNonContinuousAlphabet_CorrectAlleleMarkerEndBoundary) {
    auto prg_raw = encode_prg("7g8c8g9t10a10");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = get_allele_marker_sa_interval(8, prg_info);
    SA_Interval expected = {7, 8};
    EXPECT_EQ(result, expected);
}


/*
PRG: GCGCT5C6G6T6AGTCCT
i	BWT	SA	text_suffix
0	T	18
1	6	12	A G T C C T
2	T	15	C C T
3	G	1	C G C T 5 C 6 G 6 T 6 A G T C C T
4	C	16	C T
5	G	3	C T 5 C 6 G 6 T 6 A G T C C T
6	5	6	C 6 G 6 T 6 A G T C C T
7	0	0	G C G C T 5 C 6 G 6 T 6 A G T C C T
8	C	2	G C T 5 C 6 G 6 T 6 A G T C C T
9	A	13	G T C C T
10	6	8	G 6 T 6 A G T C C T
11	C	17	T
12	G	14	T C C T
13	C	4	T 5 C 6 G 6 T 6 A G T C C T
14	6	10	T 6 A G T C C T
15	T	5	5 C 6 G 6 T 6 A G T C C T
16	T	11	6 A G T C C T
17	C	7	6 G 6 T 6 A G T C C T
18	G	9	6 T 6 A G T C C T

 */


TEST(MarkerSearch, AtSiteEnd_GetAllMarkerChars) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
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
    std::unordered_set<uint64_t> expected = {6};
    EXPECT_EQ(result, expected);
}


TEST(MarkerSAIntervals, AtSiteExitPoint_NewSearchState_WithAllAlleles) {
    auto prg_raw = encode_prg("gcgct5c6g6t6Agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 1}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);
    EXPECT_EQ(markers_search_states.size(), 1);

    const auto first = markers_search_states.front();
    auto sa_interval = first.sa_interval;
    SA_Interval expected = SA_Interval{16,18};
    EXPECT_EQ(sa_interval, expected);
}


TEST(VariantLocus_Path, AtSiteExitPoint_VariantPathOfAllAlleles) {
    auto prg_raw = encode_prg("gcgct5c6g6t6Agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    // first char: a
    SearchState initial_search_state = {
            SA_Interval {1, 1}
    };
    const auto &markers_search_states = process_markers_search_state(initial_search_state,
                                                                     prg_info);

    std::vector<VariantLocus> result;
    for (const auto &search_state: markers_search_states)
        result.push_back(search_state.traversing_path.front());

    // We expect the following: one searchstate has two alleles, and thus the allele part is unspecified still
    // The other is a singleton SearchState corresponding to the end of the site, which is alleles #3.
    std::vector<VariantLocus> expected = {
            {5, ALLELE_UNKNOWN},
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, GivenAlleleMarkerSaIndex_ReturnAlleleId) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    uint64_t allele_marker_sa_index = 18;
    auto result = get_allele_id(allele_marker_sa_index,
                                prg_info);
    auto expected = 2;
    EXPECT_EQ(result, expected);
}


TEST(ExitASite, ThirdAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
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
            SA_Interval {15, 15},
            VariantSitePath {VariantLocus{5,3}},
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site,
    };
    EXPECT_EQ(result, expected);
}


TEST(ExitASite, SecondAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
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
            SA_Interval {15, 15},
            VariantSitePath {VariantLocus{5, 2}},
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site,
    };
    EXPECT_EQ(result, expected);
}


TEST(ExitASite, FirstAlleleSingleChar_SkipToSiteStartBoundaryMarker) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
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
            SA_Interval {15, 15},
            VariantSitePath {VariantLocus{5,1}},
            VariantSitePath {},
            SearchVariantSiteState::outside_variant_site,
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: GCGCT5C6G6T6AGTCCT
i	BWT	SA	text_suffix
0	T	18
1	6	12	A G T C C T
2	T	15	C C T
3	G	1	C G C T 5 C 6 G 6 T 6 A G T C C T
4	C	16	C T
5	G	3	C T 5 C 6 G 6 T 6 A G T C C T
6	5	6	C 6 G 6 T 6 A G T C C T
7	0	0	G C G C T 5 C 6 G 6 T 6 A G T C C T
8	C	2	G C T 5 C 6 G 6 T 6 A G T C C T
9	A	13	G T C C T
10	6	8	G 6 T 6 A G T C C T
11	C	17	T
12	G	14	T C C T
13	C	4	T 5 C 6 G 6 T 6 A G T C C T
14	6	10	T 6 A G T C C T
15	T	5	5 C 6 G 6 T 6 A G T C C T
16	T	11	6 A G T C C T
17	C	7	6 G 6 T 6 A G T C C T
18	G	9	6 T 6 A G T C C T
 */

TEST(Search, InitialStateWithPopulatedVariantSitePath_CorrectVariantSitePathInResult) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
    auto prg_info = generate_prg_info(prg_raw);
    auto pattern_char = encode_dna_base('t');

    SearchState initial_search_state = {
            SA_Interval {10,10}, //Starting at char 'c' at index 6 in prg
            VariantSitePath {},
            VariantSitePath {},
            SearchVariantSiteState::unknown,
    };
    SearchStates initial_search_states = {initial_search_state};

    auto final_search_states=process_read_char_search_states(pattern_char,initial_search_states,prg_info);

    EXPECT_EQ(final_search_states.size(), 1);
    auto search_state = final_search_states.front();
    const auto &result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 2},
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, KmerAbsentFromKmerIndex_NoSearchStatesReturned) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
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
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
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
    SA_Interval expected = {14, 14};
    EXPECT_EQ(result, expected);
}


TEST(VariantLocus_Path, GivenSearchEndingInAllele_CorrectVariantSitePath) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("tagtcc");
    Pattern kmer = encode_dna_bases("gtcc");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 3}
    };
    EXPECT_EQ(result, expected);
}


TEST(VariantLocus_Path, GivenSearchStartingInAllele_CorrectVariantSitePath) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("cgctg");
    Pattern kmer = encode_dna_bases("gctg");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(VariantLocus_Path, GivenSearchCrossingAllele_CorrectVariantSitePath) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("ctgag");
    Pattern kmer = encode_dna_bases("tgag");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: GCT5C6G6T6AG7T8C8CT
i	BWT	SA	text_suffix
0	T	19
1	6	10	A G 7 T 8 C 8 C T
2	8	17	C T
3	G	1	C T 5 C 6 G 6 T 6 A G 7 T 8 C 8 C T
4	5	4	C 6 G 6 T 6 A G 7 T 8 C 8 C T
5	8	15	C 8 C T
6	0	0	G C T 5 C 6 G 6 T 6 A G 7 T 8 C 8 C T
7	6	6	G 6 T 6 A G 7 T 8 C 8 C T
8	A	11	G 7 T 8 C 8 C T
9	C	18	T
10	C	2	T 5 C 6 G 6 T 6 A G 7 T 8 C 8 C T
11	6	8	T 6 A G 7 T 8 C 8 C T
12	7	13	T 8 C 8 C T
13	T	3	5 C 6 G 6 T 6 A G 7 T 8 C 8 C T
14	T	9	6 A G 7 T 8 C 8 C T
15	C	5	6 G 6 T 6 A G 7 T 8 C 8 C T
16	G	7	6 T 6 A G 7 T 8 C 8 C T
17	G	12	7 T 8 C 8 C T
18	C	16	8 C T
19	T	14	8 C 8 C T
*/


TEST(VariantLocus_Path, ReadCrossingTwoAlleles) {
    auto prg_raw = encode_prg("gct5c6g6t6ag7t8c8ct");
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("tct");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("cagtct");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {7, 1},
            VariantLocus {5, 1},
    };
    EXPECT_EQ(result, expected);
}


TEST(VarSiteBSearch, StartWithinAllele_MapToOtherAllele) {
    auto prg_raw = encode_prg("gct5c6g6t6ag7GAG8c8ct");
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("gag");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("caggag");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {7, 1},
            VariantLocus {5, 1},
    };
    EXPECT_EQ(result, expected);
}


TEST(VarSiteBSearch, KmerImmediatelyAfterVariantSite) {
    auto prg_raw = encode_prg("gct5c6g6t6ag7t8c8cta");
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("cta");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("gccta");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {7, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(VarSiteBSearch, KmerCrossesVariantSite) {
    auto prg_raw = encode_prg("gct5c6g6t6ag7t8c8cta");
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("gccta");
    Patterns kmers = {kmer};
    auto kmer_size = 5;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("agccta");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {7, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(EndInLocus, SearchStarts_andEnds_withinLoci) {
    auto prg_raw = encode_prg("gct5c6g6T6AG7T8c8cta");
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    auto kmer_size = 3;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("tagt");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            //The Search_State ended with ALLELE_UNKNOWN, but then we specified the allele due to completely mapping the read
            VariantLocus {7, 1},
            VariantLocus {5, 3},
    };
    EXPECT_EQ(result, expected);

    EXPECT_EQ(search_state.sa_interval.second - search_state.sa_interval.first + 1, 1);
}


/*
 * A case where we end the read mapping inside several alleles of the same site.
 * We test expected behaviour along the way from kmer indexing to read mapping alleles concurrently
 * to allele ID specification post mapping.
 */
TEST(EndInLocus, SearchEnds_AtConcurrentAlleles) {
    auto prg_raw = encode_prg("gct5gC6aC6C6t6Cg");
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
    int_Base pattern_char = 2;
    search_states = process_read_char_search_states(pattern_char,
                                                        search_states,
                                                        prg_info);

    // CONCURRENT ALLELE QUERYING
    // We expect three occurrences of 'CC' at this stage, in a single SA interval - because
    // the allele markers sort together in the SA. The allele IDs should be unspecified.
    EXPECT_EQ(search_states.size(), 1);
    EXPECT_EQ(search_states.front().traversing_path.back().second, ALLELE_UNKNOWN);

    // ALLELE ID SPECIFICATION
    // This function gets called when we have finished mapping our read and we have unknown allele ids left.
    gram::set_allele_ids(search_states, prg_info);
    EXPECT_EQ(search_states.size(), 3);

    for (const auto search_state: search_states){
        SA_Interval sa = search_states.front().sa_interval;
        EXPECT_EQ(sa.second - sa.first + 1, 1);
    }
}


TEST(VarSiteBSearch, ReadCrossesTwoVarSites) {
    auto prg_raw = encode_prg("gct5c6g6T6AG7T8c8cta");
    auto prg_info = generate_prg_info(prg_raw);

    Pattern kmer = encode_dna_bases("tagt");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto read = encode_dna_bases("cttagt");

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    EXPECT_EQ(search_states.size(), 1);

    const auto &search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {7, 1},
            VariantLocus {5, 3},
    };
    EXPECT_EQ(result, expected);
}


TEST(StartEndInLocus, OneMappingEncapsulatedByAllele) {
    auto prg_raw = encode_prg("t5c6gCTTAGT6aa");
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

    VariantLocus cov = {5, 2};
    EXPECT_EQ(search_state.traversed_path.front(), cov);
}


TEST(StartEndInLocus, TwoMappingsEncapsulatedByAllele_StateIsWithinVariantSite) {
    auto prg_raw = encode_prg("t5c6gcttagtacgcttagt6aa");
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
PRG: AC5T6CAGTAGTC6TA
i	BWT	SA	text_suffix
0	A	16
1	T	15	A
2	0	0	A C 5 T 6 C A G T A G T C 6 T A
3	C	6	A G T A G T C 6 T A
4	T	9	A G T C 6 T A
5	6	5	C A G T A G T C 6 T A
6	A	1	C 5 T 6 C A G T A G T C 6 T A
7	T	12	C 6 T A
8	A	7	G T A G T C 6 T A
9	A	10	G T C 6 T A
10	6	14	T A
11	G	8	T A G T C 6 T A
12	G	11	T C 6 T A
13	5	3	T 6 C A G T A G T C 6 T A
14	C	2	5 T 6 C A G T A G T C 6 T A
15	T	4	6 C A G T A G T C 6 T A
16	C	13	6 T A
*/
TEST(HandleAlleleEncapsulatedStates, AlleleEncapsulatedStateMissingPath_CorrectPathSet) {
    auto prg_raw = encode_prg("ac5t6cagtagtc6ta");
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
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(HandleAlleleEncapsulatedStates, AlleleEncapsulatedState_NoChange) {
    auto prg_raw = encode_prg("ac5t6cagtagtc6ta");
    auto prg_info = generate_prg_info(prg_raw);
    SearchStates search_states = {
            SearchState {
                    SA_Interval {8, 8},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    VariantSitePath {},
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
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(HandleAlleleEncapsulatedStates, SaIntervalGreaterThanOneAlleleEncapsulated_CorrectPathSet) {
    auto prg_raw = encode_prg("ac5t6cagtagtc6ta");
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
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: GCGCT5C6G6T6AGTCCT
i	BWT	SA	text_suffix
0	T	18
1	6	12	A G T C C T
2	T	15	C C T
3	G	1	C G C T 5 C 6 G 6 T 6 A G T C C T
4	C	16	C T
5	G	3	C T 5 C 6 G 6 T 6 A G T C C T
6	5	6	C 6 G 6 T 6 A G T C C T
7	0	0	G C G C T 5 C 6 G 6 T 6 A G T C C T
8	C	2	G C T 5 C 6 G 6 T 6 A G T C C T
9	A	13	G T C C T
10	6	8	G 6 T 6 A G T C C T
11	C	17	T
12	G	14	T C C T
13	C	4	T 5 C 6 G 6 T 6 A G T C C T
14	6	10	T 6 A G T C C T
15	T	5	5 C 6 G 6 T 6 A G T C C T
16	T	11	6 A G T C C T
17	C	7	6 G 6 T 6 A G T C C T
18	G	9	6 T 6 A G T C C T
*/


TEST(HandleAlleleEncapsulatedStates, OutsideSite_NoPathSet) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
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
                    VariantSitePath {},
                    SearchVariantSiteState::outside_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}

/*
PRG: CAGTAA5T6CAGTAGGC6TA
i	BWT	SA	text_suffix
0	A	20
1	T	19	A
2	T	4	A A 5 T 6 C A G T A G G C 6 T A
3	T	13	A G G C 6 T A
4	C	1	A G T A A 5 T 6 C A G T A G G C 6 T A
5	C	10	A G T A G G C 6 T A
6	A	5	A 5 T 6 C A G T A G G C 6 T A
7	0	0	C A G T A A 5 T 6 C A G T A G G C 6 T A
8	6	9	C A G T A G G C 6 T A
9	G	16	C 6 T A
10	G	15	G C 6 T A
11	A	14	G G C 6 T A
12	A	2	G T A A 5 T 6 C A G T A G G C 6 T A
13	A	11	G T A G G C 6 T A
14	6	18	T A
15	G	3	T A A 5 T 6 C A G T A G G C 6 T A
16	G	12	T A G G C 6 T A
17	5	7	T 6 C A G T A G G C 6 T A
18	A	6	5 T 6 C A G T A G G C 6 T A
19	T	8	6 C A G T A G G C 6 T A
20	C	17	6 T A
*/

TEST(HandleAlleleEncapsulatedState, ReadAlleleEncapsulatedAndOutsideSite_SplitIntoTwoSearchStates) {
    auto prg_raw = encode_prg("Cagtaa5t6Cagtaggc6ta");
    auto prg_info = generate_prg_info(prg_raw);

    SearchState search_state = {
            SA_Interval {7, 8}

    };
    auto result = handle_allele_encapsulated_state(search_state, prg_info);
    SearchStates expected = {
            SearchState {
                    SA_Interval {7, 7},
                    VariantSitePath {},
                    VariantSitePath {},
                    SearchVariantSiteState::outside_variant_site
            },
            SearchState {
                    SA_Interval {8, 8},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: TCAGTT5TCAGTCAG6ATCAGTTTCAG6TA7ATCAGT8GTG8G
i	BWT	SA	text_suffix
0	G	43
1	C	9	A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
2	C	19	A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
3	C	2	A G T T 5 T C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
4	C	34	A G T 8 G T G 8 G
5	C	13	A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
6	C	25	A G 6 T A 7 A T C A G T 8 G T G 8 G
7	6	16	A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
8	7	31	A T C A G T 8 G T G 8 G
9	T	29	A 7 A T C A G T 8 G T G 8 G
10	T	8	C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
11	T	18	C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
12	T	1	C A G T T 5 T C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
13	T	33	C A G T 8 G T G 8 G
14	T	12	C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
15	T	24	C A G 6 T A 7 A T C A G T 8 G T G 8 G
16	8	42	G
17	A	10	G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
18	8	38	G T G 8 G
19	A	20	G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
20	A	3	G T T 5 T C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
21	A	35	G T 8 G T G 8 G
22	A	14	G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
23	A	26	G 6 T A 7 A T C A G T 8 G T G 8 G
24	T	40	G 8 G
25	6	28	T A 7 A T C A G T 8 G T G 8 G
26	5	7	T C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
27	A	17	T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
28	0	0	T C A G T T 5 T C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
29	A	32	T C A G T 8 G T G 8 G
30	G	11	T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
31	T	23	T C A G 6 T A 7 A T C A G T 8 G T G 8 G
32	G	39	T G 8 G
33	T	22	T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
34	G	21	T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
35	G	4	T T 5 T C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
36	T	5	T 5 T C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
37	G	36	T 8 G T G 8 G
38	T	6	5 T C A G T C A G 6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
39	G	15	6 A T C A G T T T C A G 6 T A 7 A T C A G T 8 G T G 8 G
40	G	27	6 T A 7 A T C A G T 8 G T G 8 G
41	A	30	7 A T C A G T 8 G T G 8 G
42	G	41	8 G
43	T	37	8 G T G 8 G
*/
TEST(HandleAlleleEncapsulatedState, MappingMultipleAlleleEncapsulation_CorrectSearchStates) {
    auto prg_raw = encode_prg("tcagtt5tcagtcag6atcagtttcag6ta7atcagt8gtg8g");
    auto prg_info = generate_prg_info(prg_raw);

    // All the C's
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
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            },
            SearchState {
                    SA_Interval {11, 11},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            },
            SearchState {
                    SA_Interval {12, 12},
                    VariantSitePath {},
                    VariantSitePath {},
                    SearchVariantSiteState::outside_variant_site
            },
            SearchState {
                    SA_Interval {13, 13},
                    VariantSitePath {
                            VariantLocus {7, 1}
                    },
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            },
            SearchState {
                    SA_Interval {14, 14},
                    VariantSitePath {
                            VariantLocus {5, 1}
                    },
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            },
            SearchState {
                    SA_Interval {15, 15},
                    VariantSitePath {
                            VariantLocus {5, 2}
                    },
                    VariantSitePath {},
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(Search, ReadLeadsToPrgEdge_NoSearchStatesFound) {
    auto prg_raw = encode_prg("gcgct5c6g6t5agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("agcgc");
    Pattern kmer = encode_dna_bases("gcgc");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    ASSERT_TRUE(search_states.empty());
}
