/** @file
 * You need to distinguish tests where:
 *      - The kmer/read ends inside a variant site. Then the `traversing_path` contains latest entered site
 *      - The kmer/read ends outside a variant site. Then the `traversed_path` contains latest entered site
 */
#include <cctype>

#include "gtest/gtest.h"

#include "src_common/common.hpp"
#include "build/kmer_index/build.hpp"
#include "build/kmer_index/load.hpp"
#include "build/kmer_index/dump.hpp"


using namespace gram;


TEST(GenerateKmerIndex, GivenDnaString_DnaBasesEncodedCorrectly) {
    const auto dna_str = "AAACCCGGGTTTACGT";
    const auto result = encode_dna_bases(dna_str);
    const Sequence expected = {
            1, 1, 1,
            2, 2, 2,
            3, 3, 3,
            4, 4, 4,
            1, 2, 3, 4,
    };
    EXPECT_EQ(result, expected);
}


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


TEST(IndexKmers, KmerCrossesSecondAllele_CorrectVariantSitePath) {
    const auto prg_raw = encode_prg("acA5g6T6GCTc");
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("atgct");
    const int kmer_size = 5;
    Sequences kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);
    auto search_states = kmer_index[kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversed_path;

    VariantSitePath expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerDoesNotCrossSite_CorrectSaInterval) {
    const auto prg_raw = encode_prg("aca5g6t6gctc");
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("gctc");
    const int kmer_size = 4;
    Sequences kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);
    auto search_states = kmer_index[kmer];
    auto search_state = search_states.front();
    auto result = search_state.sa_interval;

    SA_Interval expected = {6, 6};
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerDoesNotCrossSite_CorrectVariantSitePath) {
    const auto prg_raw = encode_prg("aca5g6t6gctc");
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("gctc");
    const int kmer_size = 4;
    Sequences kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);
    auto search_states = kmer_index[kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversed_path;

    VariantSitePath expected = {};
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerCrossesFirstAllele_VariantRegionRecordedInSites) {
    const auto prg_raw = encode_prg("aca5g6t6gcatt");
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("aggca");
    const int kmer_size = 5;
    Sequences kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);
    auto search_states = kmer_index[kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversed_path;

    VariantSitePath expected = {
            VariantLocus {5, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, BothKmersOverlapVariantSiteAlleles_CorrectSearchResults) {
    auto prg_raw = encode_prg("aca5g6c6tatt");
    auto prg_info = generate_prg_info(prg_raw);

    auto kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto kmer_prefix_diff = encode_dna_bases("ac");
    Sequences kmers = {
            first_full_kmer,
            kmer_prefix_diff
    };
    auto second_full_kmer = encode_dna_bases("actat");

    auto result = index_kmers(kmers, kmer_size, prg_info);

    KmerIndex expected = {
            {first_full_kmer,
                    SearchStates {
                            SearchState {
                                    SA_Interval {3, 3},
                                    VariantSitePath {
                                            VariantLocus {5, 1}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            },
            {second_full_kmer,
                    SearchStates {
                            SearchState {
                                    SA_Interval {3, 3},
                                    VariantSitePath {
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerNotFoundInPrg_KmerAbsentFromKmerIndex) {
    auto prg_raw = encode_prg("aca5g6c6tatt");
    auto prg_info = generate_prg_info(prg_raw);

    auto kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("attat");
    auto kmer_prefix_diff = encode_dna_bases("ac");
    Sequences kmers = {
            first_full_kmer,
            kmer_prefix_diff
    };
    auto second_full_kmer = encode_dna_bases("actat");

    auto result = index_kmers(kmers, kmer_size, prg_info);

    KmerIndex expected = {
            {second_full_kmer,
                    SearchStates {
                            SearchState {
                                    SA_Interval {3, 3},
                                    VariantSitePath {
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, OneKmersOverlapsVariantSiteAllele_CorrectSearchResults) {
    const auto prg_raw = encode_prg("aca5g6c6tatt");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto kmer_prefix_diff = encode_dna_bases("aa");
    auto second_full_kmer = encode_dna_bases("aatat");
    Sequences kmers = {
            first_full_kmer,
            kmer_prefix_diff
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto first_search_states = kmer_index[first_full_kmer];
    auto first_search_state = first_search_states.front();
    auto first_result = first_search_state.traversed_path;
    VariantSitePath first_expected = {
            VariantLocus {5, 1}
    };
    EXPECT_EQ(first_result, first_expected);

    auto second_search_states = kmer_index[second_full_kmer];
    EXPECT_TRUE(second_search_states.empty());
}


TEST(IndexKmers, ThreeKmersOverlapSiteThreeAllele_CorrectSearchResults) {
    const auto prg_raw = encode_prg("aca5g6c6a6tatt");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto second_full_kmer = encode_dna_bases("actat");
    auto third_full_kmer = encode_dna_bases("aatat");
    Sequences kmers = {
            encode_dna_bases("agtat"),
            encode_dna_bases("ac"),
            encode_dna_bases("aa"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 1}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[second_full_kmer];
    search_state = search_states.front();
    result = search_state.traversed_path;
    expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[third_full_kmer];
    search_state = search_states.front();
    result = search_state.traversed_path;
    expected = {
            VariantLocus {5, 3}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, ThreeKmersOneMissMatch_CorrectSearchResults) {
    const auto prg_raw = encode_prg("aca5g6c6a6tatt");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto second_full_kmer = encode_dna_bases("actat");
    auto third_full_kmer = encode_dna_bases("attat");
    Sequences kmers = {
            encode_dna_bases("agtat"),
            encode_dna_bases("ac"),
            encode_dna_bases("at"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 1}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[second_full_kmer];
    search_state = search_states.front();
    result = search_state.traversed_path;
    expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[third_full_kmer];
    EXPECT_TRUE(search_states.empty());
}


TEST(IndexKmers, OneKmerStartsAtAllele_SiteFound) {
    const auto prg_raw = encode_prg("aca5g6c6a6tatt");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("gtat");
    Sequences kmers = {
            encode_dna_bases("gtat"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversing_path;
    VariantSitePath expected = {
            VariantLocus {5, ALLELE_UNKNOWN}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerFromAlleleCenter_KmerEntryFoundNoVariantSitePath) {
    const auto prg_raw = encode_prg("gct5cccc6g6t6ag");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 3;
    auto first_full_kmer = encode_dna_bases("ccc");
    Sequences kmers = {
            encode_dna_bases("ccc"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto found = kmer_index.find(first_full_kmer) != kmer_index.end();
    EXPECT_TRUE(found);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {};
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, TwoKmersStartAtAllele_SitesFound) {
    const auto prg_raw = encode_prg("aca5g6c6a6tatt");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("gtat");
    auto second_full_kmer = encode_dna_bases("ctat");
    // Only writing 'c' as second kmer in list
    // below means we will index 'ctat' because of prefix diffing.
    Sequences kmers = {
            encode_dna_bases("gtat"),
            encode_dna_bases("c"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversing_path;
    VariantSitePath expected = {
            VariantLocus {5, ALLELE_UNKNOWN}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[second_full_kmer];
    search_state = search_states.front();
    result = search_state.traversing_path;
    expected = {
            VariantLocus {5, ALLELE_UNKNOWN}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerEndingInAllele_SingleSiteFound) {
    const auto prg_raw = encode_prg("aca5g6c6t");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("acag");
    Sequences kmers = {
            encode_dna_bases("acag"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, TwoKmersEndingInAlleles_TwoSingleSitesFound) {
    const auto prg_raw = encode_prg("aca5g6c6t");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("acag");
    auto second_full_kmer = encode_dna_bases("acac");
    Sequences kmers = {
            encode_dna_bases("acag"),
            encode_dna_bases("acac"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 1}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[second_full_kmer];
    search_state = search_states.back();
    result = search_state.traversed_path;
    expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerStartingInSiteAndEndInAnotherSite_CorrectVariantSitePath) {
    const auto prg_raw = encode_prg("aca5g6C6TT7A8c8gg");
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("ctta");
    Sequences kmers = {
            encode_dna_bases("ctta"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = std::make_pair(search_state.traversed_path, search_state.traversing_path);
    auto expected = std::make_pair(
            std::vector{VariantLocus {7, 1} },
            std::vector{ VariantLocus {5, ALLELE_UNKNOWN} }
            );
    EXPECT_EQ(result, expected);
}


/*
PRG: TTT5TA6T6ACG
i	BWT	SA	text_suffix
0	G	12
1	6	9	A C G
2	T	5	A 6 T 6 A C G
3	A	10	C G
4	C	11	G
5	5	4	T A 6 T 6 A C G
6	0	0	T T T 5 T A 6 T 6 A C G
7	T	1	T T 5 T A 6 T 6 A C G
8	T	2	T 5 T A 6 T 6 A C G
9	6	7	T 6 A C G
10	T	3	5 T A 6 T 6 A C G
11	T	8	6 A C G
12	A	6	6 T 6 A C G
*/

TEST(IndexKmers, TwoSearchStatesIdenticalSaIntervals_DifferentVariantSitePaths) {
    auto prg_raw = encode_prg("ttt5ta6t6acg");
    auto prg_info = generate_prg_info(prg_raw);

    auto kmer_size = 4;
    auto kmer = encode_dna_bases("tttt");
    Sequences kmers = {kmer};

    auto result = index_kmers(kmers, kmer_size, prg_info);
    // Note for the expectation: the markers get processed in reverse SA index ordering
    KmerIndex expected = {
            {kmer,
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, 1}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, GivenPrgWithLongNonVariantTail_KmerEndingAtTailExtracted) {
    //               |                               |
    auto prg_raw = encode_prg("atggaacggct25cg26cc26tg26tc26cg27g28a28tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat");
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 15;
    parameters.max_read_size = 20;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {4, 3, 3, 1, 1, 2, 3, 3, 2, 4, 2, 3, 2, 3, 3};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, GivenPrgWithLongNonVariantTail_KmerStartingAtLeftMostAlleleCharExtracted) {
    //                                                  |                          |
    auto prg_raw = encode_prg("atggaacggct25cg26cc26tg26tc26cg27g28a28tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat");
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 15;
    parameters.max_read_size = 20;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {1, 4, 2, 2, 2, 2, 3, 1, 2, 3, 1, 4, 4, 2, 2};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, GivenPrgWithLongNonVariantTail_KmerImmediatelyAfterSiteExtracted) {
    //                                                     |                        |
    auto prg_raw = encode_prg("atggaacggct25cg26cc26tg26tc26cg27g28a28tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat");
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 15;
    parameters.max_read_size = 20;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {4, 2, 2, 2, 2, 3, 1, 2, 3, 1, 4, 4, 2, 2, 2};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, KmerStartsOneBaseBeyondRangeEdge_KmerNotExtracted) {
    //                                                                     |             |
    auto prg_raw = encode_prg("atggaacggct25cg26cc26tg26tc26cg27g28a28tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat");
    //                                                                    ^region end
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 15;
    parameters.max_read_size = 20;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {3, 1, 2, 3, 1, 4, 4, 2, 2, 2, 2, 3, 1, 2, 3};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_FALSE(found);
}


TEST(IndexKmers, KmerStartsAtRangeEdge_KmerExtracted) {
    //                                                                     |             |
    auto prg_raw = encode_prg("atggaacggct25cg26cc26tg26tc26cg27g28a28tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat");
    //                                                                     ^region end
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 15;
    parameters.max_read_size = 21;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {3, 1, 2, 3, 1, 4, 4, 2, 2, 2, 2, 3, 1, 2, 3};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, KmerWithinMaxReadSizeRegionNoSiteOverlap_KmerFound) {
    //                 last site overlapping kmer end: |
    auto prg_raw = encode_prg("t25cg26cc26tg26tc26ctcacagacgattctcctgac");
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 18;
    parameters.max_read_size = 22;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {1, 2, 1, 3, 1, 2, 3, 1, 4, 4, 2, 4, 2, 2, 4, 3, 1, 2};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, KmerEndJustOutsideMaxReadSize_KmerNotFoundInIndex) {
    //                 last site overlapping kmer end: |
    auto prg_raw = encode_prg("t25cg26cc26tg26tc26ctcacagacgattctcctgac");
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 18;
    parameters.max_read_size = 21;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {1, 2, 1, 3, 1, 2, 3, 1, 4, 4, 2, 4, 2, 2, 4, 3, 1, 2};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_FALSE(found);
}


TEST(IndexKmers, TwoSitesAndKmerWithinMaxReadSizeRegionNoSiteOverlap_KmerFound) {
    //                  last base given max read size:   |
    auto prg_raw = encode_prg("t25cg26cc26tg26tc26ct27ca28ca28gacgattctcctgac");
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 5;
    parameters.max_read_size = 8;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {2, 3, 1, 4, 4};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, TwoSitesAndKmerOutsideMaxReadSizeRegionNoSiteOverlap_KmerNotFound) {
    //                  last base given max read size:   |
    auto prg_raw = encode_prg("t25cg26cc26tg26tc26ct27ca28ca28gacgattctcctgac");
    auto prg_info = generate_prg_info(prg_raw);

    BuildParams parameters = {};
    parameters.kmers_size = 12;
    parameters.max_read_size = 7;

    auto kmer_prefix_diffs = get_all_kmer_and_compute_prefix_diffs(parameters,
                                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Sequence target_kmer = {2, 3, 1, 4, 4, 2, 4, 2, 2, 4, 3, 1};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_FALSE(found);
}


TEST(IndexKmers, GivenTwoSerializedKmers_CorrectlyExtrctedKmers) {
    sdsl::int_vector<3> all_kmers = {1, 2, 3, 4, 1, 2, 1, 2};
    const uint32_t kmer_size = 4;
    uint64_t kmer_start_index = 0;

    std::vector<Sequence> result = {};
    auto kmer = deserialize_next_kmer(kmer_start_index, all_kmers, kmer_size);
    result.push_back(kmer);

    kmer_start_index += kmer_size;
    kmer = deserialize_next_kmer(kmer_start_index, all_kmers, kmer_size);
    result.push_back(kmer);

    std::vector<Sequence> expected = {
            {1, 2, 3, 4},
            {1, 2, 1, 2},
    };
    EXPECT_EQ(result, expected);
}

