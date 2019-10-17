#include "gtest/gtest.h"
#include "src_common/generate_prg.hpp"
#include "quasimap/search/encapsulated_search.hpp"

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
            SearchState{
                    SA_Interval{8, 8}
            }
    };
    auto result = handle_allele_encapsulated_states(search_states, prg_info);
    SearchStates expected = {
            SearchState{
                    SA_Interval{8, 8},
                    VariantSitePath{
                            VariantLocus{5, 2}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(HandleAlleleEncapsulatedStates, AlleleEncapsulatedState_NoChange) {
    auto prg_raw = encode_prg("ac5t6cagtagtc6ta");
    auto prg_info = generate_prg_info(prg_raw);
    SearchStates search_states = {
            SearchState{
                    SA_Interval{8, 8},
                    VariantSitePath{
                            VariantLocus{5, 2}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            }
    };
    auto result = handle_allele_encapsulated_states(search_states, prg_info);
    SearchStates expected = {
            SearchState{
                    SA_Interval{8, 8},
                    VariantSitePath{
                            VariantLocus{5, 2}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(HandleAlleleEncapsulatedStates, SaIntervalGreaterThanOneAlleleEncapsulated_CorrectPathSet) {
    auto prg_raw = encode_prg("ac5t6cagtagtc6ta");
    auto prg_info = generate_prg_info(prg_raw);
    SearchStates search_states = {
            SearchState{
                    SA_Interval{3, 4} // interval of all 'ag'
            }
    };
    auto result = handle_allele_encapsulated_states(search_states, prg_info);
    // Both mappings are encapsulated in same site & allele, so have a size 2 SA interval.
    SearchStates expected = {
            SearchState{
                    SA_Interval{3, 4},
                    VariantSitePath{
                            VariantLocus{5, 2}
                    },
                    VariantSitePath{},
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
            SearchState{
                    SA_Interval{7, 7}
            }
    };
    auto result = handle_allele_encapsulated_states(search_states, prg_info);
    SearchStates expected = {
            SearchState{
                    SA_Interval{7, 7},
                    VariantSitePath{},
                    VariantSitePath{},
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
            SA_Interval{7, 8}

    };
    auto result = handle_allele_encapsulated_state(search_state, prg_info);
    SearchStates expected = {
            SearchState{
                    SA_Interval{7, 7},
                    VariantSitePath{},
                    VariantSitePath{},
                    SearchVariantSiteState::outside_variant_site
            },
            SearchState{
                    SA_Interval{8, 8},
                    VariantSitePath{
                            VariantLocus{5, 2}
                    },
                    VariantSitePath{},
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
            SA_Interval{10, 15}

    };
    auto result = handle_allele_encapsulated_state(search_state, prg_info);
    SearchStates expected = {
            SearchState{
                    SA_Interval{10, 10},
                    VariantSitePath{
                            VariantLocus{5, 1}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            },
            SearchState{
                    SA_Interval{11, 11},
                    VariantSitePath{
                            VariantLocus{5, 2}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            },
            SearchState{
                    SA_Interval{12, 12},
                    VariantSitePath{},
                    VariantSitePath{},
                    SearchVariantSiteState::outside_variant_site
            },
            SearchState{
                    SA_Interval{13, 13},
                    VariantSitePath{
                            VariantLocus{7, 1}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            },
            SearchState{
                    SA_Interval{14, 14},
                    VariantSitePath{
                            VariantLocus{5, 1}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            },
            SearchState{
                    SA_Interval{15, 15},
                    VariantSitePath{
                            VariantLocus{5, 2}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            }
    };
    EXPECT_EQ(result, expected);
}


