/**
 * @file
 * Unit tests for regular backward searching.
 *
 * Test suites:
 *  - NoVarPrg: backward searching on Prg with no variant sites.
 *  - VarPrg: backward searching on non-nested Prg with variant sites.
 *  - HandleAlleleEncapsulatedState: handle read fully mapping inside a site
 */
#include "gtest/gtest.h"
#include "src_common/submod_resources.hpp"
#include "prg/prg_info.hpp"
#include "genotype/quasimap/search/BWT_search.hpp"
#include "genotype/quasimap/quasimap.hpp"
#include "build/kmer_index/build.hpp"


/*
PRG: ACA5G6T6GCTC
i	BWT	SA	text_suffix
0	C	12
1	0	0	A C A 5 G 6 T 6 G C T C
2	C	2	A 5 G 6 T 6 G C T C
3	T	11	C
4	A	1	C A 5 G 6 T 6 G C T C
5	G	9	C T C
6	6	8	G C T C
7	5	4	G 6 T 6 G C T C
8	C	10	T C
9	6	6	T 6 G C T C
10	A	3	5 G 6 T 6 G C T C
11	T	7	6 G C T C
12	G	5	6 T 6 G C T C
*/

TEST(BWT_DNA_masks, rankQueries){
    const auto prg_raw = encode_prg("aca5g6t6gctc");
    auto prg_info = generate_prg_info(prg_raw);
    // The interval is all suffixes starting with 'T'
    int sa_start = 8;
    int sa_end = 9;
    // How many 'C' up to and excluding sa_start?
    EXPECT_EQ(gram::dna_bwt_rank(sa_start, 2, prg_info), 2);
    // How many 'C' up to and including sa_end?
    EXPECT_EQ(gram::dna_bwt_rank(sa_end, 2, prg_info), 3);
}

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

TEST(noVarPrg, SingleChar_CorrectSaIntervalReturned) {
    auto prg_raw = encode_prg("gcgctggagtgctgt");
    auto prg_info = generate_prg_info(prg_raw);
    auto pattern_char = encode_dna_base('g');

    SearchState initial_search_state = {
            SA_Interval{0, prg_info.fm_index.size() - 1}
    };
    SearchStates search_states = {initial_search_state};

    auto result = search_base_backwards(pattern_char,
                                        search_states,
                                        prg_info);
    SearchStates expected = {
            SearchState{
                    SA_Interval{5, 11},
                    VariantSitePath{}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(noVarPrg, TwoConsecutiveChars_CorrectFinalSaIntervalReturned) {
    auto prg_raw = encode_prg("gcgctggagtgctgt");
    auto prg_info = generate_prg_info(prg_raw);

    SearchState initial_search_state = {
            SA_Interval{0, prg_info.fm_index.size() - 1}
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
            SearchState{
                    SA_Interval{13, 15},
                    VariantSitePath{}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(noVarPrg, SingleCharFreqOneInText_SingleSA) {
    auto prg_raw = encode_prg("gcgctggagtgctgt");
    auto prg_info = generate_prg_info(prg_raw);
    auto pattern_char = encode_dna_base('a');

    SearchState initial_search_state = {
            SA_Interval{0, prg_info.fm_index.size() - 1}
    };
    SearchStates search_states = {initial_search_state};

    auto result = search_base_backwards(pattern_char,
                                        search_states,
                                        prg_info);
    SearchStates expected = {
            SearchState{
                    SA_Interval{1, 1},
                    VariantSitePath{}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(noVarPrg, TwoConsecutiveChars_SingleSaIntervalEntry) {
    auto prg_raw = encode_prg("gcgctggagtgctgt");
    auto prg_info = generate_prg_info(prg_raw);

    SearchState initial_search_state = {
            SA_Interval{0, prg_info.fm_index.size() - 1}
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


TEST(noVarPrg, TwoConsecutiveCharsNoValidSaInterval_NoSearchStatesReturned) {
    auto prg_raw = encode_prg("gcgctggagtgctgt");
    auto prg_info = generate_prg_info(prg_raw);

    SearchState initial_search_state = {
            SA_Interval{0, prg_info.fm_index.size() - 1}
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

TEST(VarPrg, oneBaseExtensionGC_CorrectSaInterval) {
    // Looking for 'GC' here
    auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    Marker next_char = 3;
    SA_Index next_char_first_sa_index = 8; // Where the first 'G' lies
    SA_Interval current_sa_interval = {3, 7}; // start at 'C'

    auto result = base_next_sa_interval(next_char,
                                        next_char_first_sa_index,
                                        current_sa_interval,
                                        prg_info);
    SA_Interval expected = {8, 9};
    EXPECT_EQ(result, expected);
}


TEST(VarPrg, oneBaseExtensionAG_CorrectSaInterval) {
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


TEST(VarPrg, ReadLeadsToPrgEdge_NoSearchStatesFound) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    auto read = encode_dna_bases("agcgc");
    Sequence kmer = encode_dna_bases("gcgc");
    Sequences kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    ASSERT_TRUE(search_states.empty());
}
