#include <cctype>

#include "gtest/gtest.h"

#include "kmer_index/kmer_index.hpp"
#include "../test_utils.hpp"


TEST(GenerateKmerIndex, GivenDataForSingleKmerIndexEntry_CorrectRowDumpGenerated) {
    VariantSitePath first_path = {
            VariantSite {5, 9},
            VariantSite {7, 19},
            VariantSite {9, 1},
    };
    VariantSitePath second_path = {
            VariantSite {9, 29},
            VariantSite {11, 39},
    };

    SearchStates search_states = {
            SearchState {
                    SA_Interval {123, 456},
                    first_path
            },
            SearchState {
                    SA_Interval {789, 424},
                    second_path
            }
    };

    const auto result = dump_kmer_index_entry(search_states);
    const auto expected = "123 456 789 424|5 9 7 19 9 1|9 29 11 39|";
    EXPECT_EQ(result, expected);
}


TEST(GenerateKmerIndex, TwoSearchStateOneVairiantPath_CorrectKmerIndexEntryDump) {
    VariantSitePath first_path = {
            VariantSite {5, 9},
            VariantSite {7, 19},
            VariantSite {9, 1},
    };
    VariantSitePath second_path = {
            VariantSite {9, 29},
            VariantSite {11, 39},
    };

    SearchStates search_states = {
            SearchState {
                    SA_Interval {123, 456},
                    first_path
            },
            SearchState {
                    SA_Interval {789, 424},
            }
    };

    const auto result = dump_kmer_index_entry(search_states);
    const auto expected = "123 456 789 424|5 9 7 19 9 1||";
    EXPECT_EQ(result, expected);
}


TEST(GenerateKmerIndex, GivenVariantSitePaths_DumpVariantSitesPathsCorrectly) {
    VariantSitePath first_path = {
            VariantSite {5, 9},
            VariantSite {7, 19},
            VariantSite {9, 1},
    };
    VariantSitePath second_path = {
            VariantSite {9, 29},
            VariantSite {11, 39},
    };

    SearchStates search_states = {
            SearchState {
                    SA_Interval {},
                    first_path
            },
            SearchState {
                    SA_Interval {},
                    second_path
            }
    };
    const auto result = dump_variant_site_paths(search_states);
    const auto expected = "5 9 7 19 9 1|9 29 11 39|";
    EXPECT_EQ(result, expected);
}


TEST(GenerateKmerIndex, GivenSaIntervals_DumpSaIntervalsStringCorrectly) {
    SearchStates search_states = {
            SearchState {
                    SA_Interval {1, 2}
            },
            SearchState {
                    SA_Interval {3, 4}
            }
    };

    auto result = dump_sa_intervals(search_states);
    const auto expected = "1 2 3 4";
    EXPECT_EQ(result, expected);
}


TEST(GenerateKmerIndex, GivenKmer_DumpKmerStringCorrectly) {
    const Pattern kmer = {1, 2, 3, 4};
    const auto result = dump_kmer(kmer);
    const auto expected = "1 2 3 4";
    EXPECT_EQ(result, expected);
}


TEST(GenerateKmerIndex, GivenDnaString_DnaBasesEncodedCorrectly) {
    const auto dna_str = "AAACCCGGGTTTACGT";
    const auto result = encode_dna_bases(dna_str);
    const Pattern expected = {
            1, 1, 1,
            2, 2, 2,
            3, 3, 3,
            4, 4, 4,
            1, 2, 3, 4,
    };
    EXPECT_EQ(result, expected);
}


TEST(ParseKmerIndex, GivenKmerIndexEntryStr_CorrectlyParsed) {
    const auto entry = "123 456 789 424|5 9 7 19 9 1|9 29 11 39|";
    auto result = parse_kmer_index_entry(entry);
    SearchStates expected = {
            SearchState {
                    SA_Interval {123, 456},
                    VariantSitePath {
                            VariantSite {5, 9},
                            VariantSite {7, 19},
                            VariantSite {9, 1}
                    }
            },

            SearchState {
                    SA_Interval {789, 424},
                    VariantSitePath {
                            VariantSite {9, 29},
                            VariantSite {11, 39}
                    }
            }

    };
    EXPECT_EQ(result, expected);
}


TEST(ParseKmerIndex, IndexEntryTwoSearchStatesOneVariantSitePath_ParsedCorrectly) {
    KmerIndex kmer_index;
    const auto entry = "123 456 789 424||9 29 11 39|";

    auto result = parse_kmer_index_entry(entry);
    SearchStates expected = {
            SearchState {
                    SA_Interval {123, 456},
                    VariantSitePath {}
            },

            SearchState {
                    SA_Interval {789, 424},
                    VariantSitePath {
                            VariantSite {9, 29},
                            VariantSite {11, 39}
                    }
            }

    };
    EXPECT_EQ(result, expected);
}


TEST(ParseKmerIndex, GivenSaIntervalsString_CorrectlyParsed) {
    const auto full_sa_intervals_str = "352511 352512 352648 352649 2 3";
    const auto result = parse_sa_intervals(full_sa_intervals_str);

    std::vector<SA_Interval> expected{
            SA_Interval {352511, 352512},
            SA_Interval {352648, 352649},
            SA_Interval {2, 3},
    };
    EXPECT_EQ(result, expected);
}


TEST(ParseKmerIndex, GivenTwoSites_CorrectSiteStructGenerated) {
    const auto kmer_index_entry = "5 9 7 19";
    const auto &result = parse_variant_site_path(kmer_index_entry);
    VariantSitePath expected = {
            VariantSite {5, 9},
            VariantSite {7, 19},
    };
    EXPECT_EQ(result, expected);
}


TEST(ParseKmerIndex, GivenSitesTrailingAt_TrailingAtIgnored) {
    const auto kmer_index_entry = "5 9 7 19";
    const auto &result = parse_variant_site_path(kmer_index_entry);
    VariantSitePath expected = {
            VariantSite {5, 9},
            VariantSite {7, 19},
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: aca5g6t5gctc
i	F	BWT	text	SA	suffix
0	0	2	1	    12	  0
1	1	0	2	    0	  1 2 1 5 3 6 4 5 3 2 4 2 0
2	1	2	1	    2	  1 5 3 6 4 5 3 2 4 2 0
3	2	4	5	    11	  2 0
4	2	1	3	    1	  2 1 5 3 6 4 5 3 2 4 2 0
5	2	3	6	    9	  2 4 2 0
6	3	5	4	    8	  3 2 4 2 0
7	3	5	5	    4	  3 6 4 5 3 2 4 2 0
8	4	2	3	    10	  4 2 0
9	4	6	2	    6	  4 5 3 2 4 2 0
10	5	4	4	    7	  5 3 2 4 2 0
11	5	1	2	    3	  5 3 6 4 5 3 2 4 2 0
12	6	3	0	    5	  6 4 5 3 2 4 2 0
 */


TEST(IndexKmers, KmerCrossesSecondAllele_CorrectVariantSitePath) {
    const std::string prg_raw = "aca5g6t5gctc";
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("atgct");
    const int kmer_size = 5;
    Patterns kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);
    auto search_states = kmer_index[kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;

    VariantSitePath expected = {
            VariantSite {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerCrossesFirstAllele_VariantRegionRecordedInSites) {
    const std::string prg_raw = "aca5g6t5gcatt";
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("aggca");
    const int kmer_size = 5;
    Patterns kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);
    auto search_states = kmer_index[kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;

    VariantSitePath expected = {
            VariantSite {5, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, BothKmersOverlapVariantSiteAlleles_CorrectSearchResults) {
    auto prg_raw = "aca5g6c5tatt";
    auto prg_info = generate_prg_info(prg_raw);

    auto kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto kmer_prefix_diff = encode_dna_bases("ac");
    Patterns kmers = {
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
                                            VariantSite {5, 1}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            },
            {second_full_kmer,
                    SearchStates {
                            SearchState {
                                    SA_Interval {3, 3},
                                    VariantSitePath {
                                            VariantSite {5, 2}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerNotFoundInPrg_KmerAbsentFromKmerIndex) {
    auto prg_raw = "aca5g6c5tatt";
    auto prg_info = generate_prg_info(prg_raw);

    auto kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("attat");
    auto kmer_prefix_diff = encode_dna_bases("ac");
    Patterns kmers = {
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
                                            VariantSite {5, 2}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, OneKmersOverlapsVariantSiteAllele_CorrectSearchResults) {
    const std::string prg_raw = "aca5g6c5tatt";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto kmer_prefix_diff = encode_dna_bases("aa");
    auto second_full_kmer = encode_dna_bases("aatat");
    Patterns kmers = {
            first_full_kmer,
            kmer_prefix_diff
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto first_search_states = kmer_index[first_full_kmer];
    auto first_search_state = first_search_states.front();
    auto first_result = first_search_state.variant_site_path;
    VariantSitePath first_expected = {
            VariantSite {5, 1}
    };
    EXPECT_EQ(first_result, first_expected);

    auto second_search_states = kmer_index[second_full_kmer];
    EXPECT_TRUE(second_search_states.empty());
}


TEST(IndexKmers, ThreeKmersOverlapSiteThreeAllele_CorrectSearchResults) {
    const std::string prg_raw = "aca5g6c6a5tatt";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto second_full_kmer = encode_dna_bases("actat");
    auto third_full_kmer = encode_dna_bases("aatat");
    Patterns kmers = {
            encode_dna_bases("agtat"),
            encode_dna_bases("ac"),
            encode_dna_bases("aa"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantSite {5, 1}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[second_full_kmer];
    search_state = search_states.front();
    result = search_state.variant_site_path;
    expected = {
            VariantSite {5, 2}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[third_full_kmer];
    search_state = search_states.front();
    result = search_state.variant_site_path;
    expected = {
            VariantSite {5, 3}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, ThreeKmersOneMissMatch_CorrectSearchResults) {
    const std::string prg_raw = "aca5g6c6a5tatt";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto second_full_kmer = encode_dna_bases("actat");
    auto third_full_kmer = encode_dna_bases("attat");
    Patterns kmers = {
            encode_dna_bases("agtat"),
            encode_dna_bases("ac"),
            encode_dna_bases("at"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantSite {5, 1}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[second_full_kmer];
    search_state = search_states.front();
    result = search_state.variant_site_path;
    expected = {
            VariantSite {5, 2}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[third_full_kmer];
    EXPECT_TRUE(search_states.empty());
}


TEST(IndexKmers, OneKmerStartsAtAllele_SiteFound) {
    const std::string prg_raw = "aca5g6c6a5tatt";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("gtat");
    Patterns kmers = {
            encode_dna_bases("gtat"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantSite {5, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerFromAlleleCenter_KmerEntryFoundNoVariantSitePath) {
    const std::string prg_raw = "gct5cccc6g6t5ag";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 3;
    auto first_full_kmer = encode_dna_bases("ccc");
    Patterns kmers = {
            encode_dna_bases("ccc"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto found = kmer_index.find(first_full_kmer) != kmer_index.end();
    EXPECT_TRUE(found);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {};
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, TwoKmersStartAtAllele_SitesFound) {
    const std::string prg_raw = "aca5g6c6a5tatt";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("gtat");
    auto second_full_kmer = encode_dna_bases("ctat");
    Patterns kmers = {
            encode_dna_bases("gtat"),
            encode_dna_bases("c"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantSite {5, 1}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[second_full_kmer];
    search_state = search_states.front();
    result = search_state.variant_site_path;
    expected = {
            VariantSite {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerEndingInAllele_SingleSiteFound) {
    const std::string prg_raw = "aca5g6c5t";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("acag");
    Patterns kmers = {
            encode_dna_bases("acag"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantSite {5, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, TwoKmersEndingInAlleles_TwoSingleSitesFound) {
    const std::string prg_raw = "aca5g6c5t";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("acag");
    auto second_full_kmer = encode_dna_bases("acac");
    Patterns kmers = {
            encode_dna_bases("acag"),
            encode_dna_bases("acac"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantSite {5, 1}
    };
    EXPECT_EQ(result, expected);

    search_states = kmer_index[second_full_kmer];
    search_state = search_states.front();
    result = search_state.variant_site_path;
    expected = {
            VariantSite {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, KmerStartingInSiteAndEndInAnotherSite_CorrectVariantSitePath) {
    const std::string prg_raw = "aca5g6c5tt7a8c7gg";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("ctta");
    Patterns kmers = {
            encode_dna_bases("ctta"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto search_states = kmer_index[first_full_kmer];
    auto search_state = search_states.front();
    auto result = search_state.variant_site_path;
    VariantSitePath expected = {
            VariantSite {5, 2},
            VariantSite {7, 1}
    };
    EXPECT_EQ(result, expected);
}


/*
PRG: ttt5ta6t5acg
i	F	BWT	text	SA	suffix
0	0	3	4	    12	0
1	1	5	4	    9	1 2 3 0
2	1	4	4	    5	1 6 4 5 1 2 3 0
3	2	1	5	    10	2 3 0
4	3	2	4	    11	3 0
5	4	5	1	    4	4 1 6 4 5 1 2 3 0
6	4	0	6	    0	4 4 4 5 4 1 6 4 5 1 2 3 0
7	4	4	4	    1	4 4 5 4 1 6 4 5 1 2 3 0
8	4	6	5	    7	4 5 1 2 3 0
9	4	4	1	    2	4 5 4 1 6 4 5 1 2 3 0
10	5	4	2	    8	5 1 2 3 0
11	5	4	3	    3	5 4 1 6 4 5 1 2 3 0
12	6	1	0	    6	6 4 5 1 2 3 0
*/
TEST(IndexKmers, TwoSearchStatesIdenticalSaIntervals_DifferentVariantSitePaths) {
    auto prg_raw = "ttt5ta6t5acg";
    auto prg_info = generate_prg_info(prg_raw);

    auto kmer_size = 4;
    auto kmer = encode_dna_bases("tttt");
    Patterns kmers = {kmer};

    auto result = index_kmers(kmers, kmer_size, prg_info);
    KmerIndex expected = {
            {kmer,
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantSite {5, 1}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantSite {5, 2}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(IndexKmers, GivenPrgWithLongNonVariantTail_KmerEndingAtTailExtracted) {
    //               |                               |
    auto prg_raw = "atggaacggct25cg26cc26tg26tc25cg27g28a27tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat";
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 15;
    parameters.max_read_size = 20;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {4, 3, 3, 1, 1, 2, 3, 3, 2, 4, 2, 3, 2, 3, 3};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, GivenPrgWithLongNonVariantTail_KmerStartingAtLeftMostAlleleCharExtracted) {
    //                                                  |                          |
    auto prg_raw = "atggaacggct25cg26cc26tg26tc25cg27g28a27tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat";
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 15;
    parameters.max_read_size = 20;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {1, 4, 2, 2, 2, 2, 3, 1, 2, 3, 1, 4, 4, 2, 2};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, GivenPrgWithLongNonVariantTail_KmerImmediatelyAfterSiteExtracted) {
    //                                                     |                        |
    auto prg_raw = "atggaacggct25cg26cc26tg26tc25cg27g28a27tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat";
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 15;
    parameters.max_read_size = 20;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {4, 2, 2, 2, 2, 3, 1, 2, 3, 1, 4, 4, 2, 2, 2};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, KmerStartsOneBaseBeyondRangeEdge_KmerNotExtracted) {
    //                                                                     |             |
    auto prg_raw = "atggaacggct25cg26cc26tg26tc25cg27g28a27tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat";
    //                                                                    ^region end
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 15;
    parameters.max_read_size = 20;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {3, 1, 2, 3, 1, 4, 4, 2, 2, 2, 2, 3, 1, 2, 3};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_FALSE(found);
}


TEST(IndexKmers, KmerStartsAtRangeEdge_KmerExtracted) {
    //                                                                     |             |
    auto prg_raw = "atggaacggct25cg26cc26tg26tc25cg27g28a27tccccgacgattccccgacgattccccgacgattccccgacgattccccgacgattccccgacgat";
    //                                                                     ^region end
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 15;
    parameters.max_read_size = 21;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {3, 1, 2, 3, 1, 4, 4, 2, 2, 2, 2, 3, 1, 2, 3};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, KmerWithinMaxReadSizeRegionNoSiteOverlap_KmerFound) {
    //                 last site overlapping kmer end: |
    auto prg_raw = "t25cg26cc26tg26tc25ctcacagacgattctcctgac";
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 18;
    parameters.max_read_size = 22;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {1, 2, 1, 3, 1, 2, 3, 1, 4, 4, 2, 4, 2, 2, 4, 3, 1, 2};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, KmerEndJustOutsideMaxReadSize_KmerNotFoundInIndex) {
    //                 last site overlapping kmer end: |
    auto prg_raw = "t25cg26cc26tg26tc25ctcacagacgattctcctgac";
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 18;
    parameters.max_read_size = 21;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {1, 2, 1, 3, 1, 2, 3, 1, 4, 4, 2, 4, 2, 2, 4, 3, 1, 2};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_FALSE(found);
}


TEST(IndexKmers, TwoSitesAndKmerWithinMaxReadSizeRegionNoSiteOverlap_KmerFound) {
    //                  last base given max read size:   |
    auto prg_raw = "t25cg26cc26tg26tc25ct27ca28ca27gacgattctcctgac";
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 5;
    parameters.max_read_size = 8;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {2, 3, 1, 4, 4};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_TRUE(found);
}


TEST(IndexKmers, TwoSitesAndKmerOutsideMaxReadSizeRegionNoSiteOverlap_KmerNotFound) {
    //                  last base given max read size:   |
    auto prg_raw = "t25cg26cc26tg26tc25ct27ca28ca27gacgattctcctgac";
    auto prg_info = generate_prg_info(prg_raw);

    Parameters parameters;
    parameters.kmers_size = 5;
    parameters.max_read_size = 7;

    auto kmer_prefix_diffs = get_kmer_prefix_diffs(parameters,
                                                   prg_info);
    auto kmer_index = index_kmers(kmer_prefix_diffs,
                                  parameters.kmers_size,
                                  prg_info);
    Pattern target_kmer = {2, 3, 1, 4, 4};
    auto found = kmer_index.find(target_kmer) != kmer_index.end();
    EXPECT_FALSE(found);
}


TEST(IndexKmers, GivenTwoSerializedKmers_CorrectlyExtrctedKmers) {
    sdsl::int_vector<3> all_kmers = {1, 2, 3, 4, 1, 2, 1, 2};
    const uint32_t kmer_size = 4;
    uint64_t kmer_start_index = 0;

    std::vector<Pattern> result = {};
    auto kmer = deserialize_next_kmer(kmer_start_index, all_kmers, kmer_size);
    result.push_back(kmer);

    kmer_start_index += kmer_size;
    kmer = deserialize_next_kmer(kmer_start_index, all_kmers, kmer_size);
    result.push_back(kmer);

    std::vector<Pattern> expected = {
            {1, 2, 3, 4},
            {1, 2, 1, 2},
    };
    EXPECT_EQ(result, expected);
}