#include <cctype>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "gtest/gtest.h"

#include "utils.hpp"
#include "prg.hpp"
#include "kmers.hpp"


TEST(GeneratePrecalc, GivenDataForSinglePrecalcEntry_CorrectDumpRowGenerated) {
    const Pattern kmer = {1, 2, 3, 4};
    const NonSiteCrossingKmers non_site_crossing_kmers = {kmer};

    VariantSitePath first_path = {
            VariantSite {5, 9},
            VariantSite {7, 19},
            VariantSite {9, 1},
    };
    VariantSitePath second_path = {
            VariantSite {9, 29},
            VariantSite {11, 39},
    };

    VariantSitePaths paths = {
            first_path,
            second_path
    };
    const KmerVariantSitePaths kmer_sites = {
            {kmer, paths},
    };

    const SA_Intervals sa_intervals = {
            SA_Interval {123, 456},
            SA_Interval {789, 424}
    };

    const auto result = dump_kmer_index_entry(kmer,
                                              sa_intervals,
                                              kmer_sites,
                                              non_site_crossing_kmers);
    const auto expected = "1 2 3 4|1|123 456 789 424||5 9 @7 19 @9 1 @|9 29 @11 39 @|";
    EXPECT_EQ(result, expected);
}


TEST(GeneratePrecalc, GivenSites_DumpSitesCorrectly) {
    VariantSitePath first_path = {
            VariantSite {5, 9},
            VariantSite {7, 19},
            VariantSite {9, 1},
    };
    VariantSitePath second_path = {
            VariantSite {9, 29},
            VariantSite {11, 39},
    };

    VariantSitePaths paths = {
            first_path,
            second_path
    };

    const Pattern kmer = {1, 2, 3, 4};
    const KmerVariantSitePaths kmer_sites = {
            {kmer, paths},
    };
    const auto result = dump_variant_site_paths(kmer, kmer_sites);
    const auto expected = "5 9 @7 19 @9 1 @|9 29 @11 39 @|";
    EXPECT_EQ(result, expected);
}


TEST(GeneratePrecalc, GivenSaIntervals_DumpSaIntervalsStringCorrectly) {
    SA_Intervals sa_intervals = {
            SA_Interval {1, 2},
            SA_Interval {3, 4}
    };
    const auto result = dump_sa_intervals(sa_intervals);
    const auto expected = "1 2 3 4";
    EXPECT_EQ(result, expected);
}


TEST(GeneratePrecalc, GivenKmer_DumpKmerStringCorrectly) {
    const Pattern kmer = {1, 2, 3, 4};
    const auto result = dump_kmer(kmer);
    const auto expected = "1 2 3 4";
    EXPECT_EQ(result, expected);
}


TEST(GeneratePrecalc, GivenDnaString_DnaBasesEncodedCorrectly) {
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


TEST(ParsePrecalc, GivenKmerIndexEntryStr_SaIntervalsParsedCorrectly) {
    KmerIndex kmer_index;
    const auto entry = "1 2 3 4|1|123 456 789 424||5 9 @7 19 @9 1 @|9 29 @11 39 @|";
    parse_kmer_index_entry(kmer_index, entry);

    const auto &result = kmer_index.sa_intervals_map;
    KmerSA_Intervals expected = {
            {Pattern {1, 2, 3, 4},
                    SA_Intervals {
                            SA_Interval {123, 456},
                            SA_Interval {789, 424}
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenKmerIndexEntryStr_VariantSitePathsCorrectlyParsed) {
    KmerIndex kmer_index;
    const auto entry = "1 2 3 4|1|123 456 789 424||5 9 @7 19 @9 1 @|9 29 @11 39 @|";
    parse_kmer_index_entry(kmer_index, entry);

    const auto &result = kmer_index.variant_site_paths_map;
    KmerVariantSitePaths expected = {
            {Pattern {1, 2, 3, 4},
                    VariantSitePaths {
                            VariantSitePath {
                                    VariantSite {5, 9},
                                    VariantSite {7, 19},
                                    VariantSite {9, 1}
                            },

                            VariantSitePath {
                                    VariantSite {9, 29},
                                    VariantSite {11, 39}
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenEncodedKmerString_CorrectlyParsed) {
    const auto encoded_kmer_str = "3 4 2 1 1 3 1 1 2";
    const auto result = parse_encoded_kmer(encoded_kmer_str);
    const Pattern expected = {3, 4, 2, 1, 1, 3, 1, 1, 2};
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenSaIntervalsString_CorrectlyParsed) {
    const auto full_sa_intervals_str = "352511 352512 352648 352649 352648 352649";
    const auto result = parse_sa_intervals(full_sa_intervals_str);

    SA_Intervals expected{
            SA_Interval {352511, 352512},
            SA_Interval {352648, 352649},
            SA_Interval {352648, 352649},
    };
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenTwoSites_CorrectSiteStructGenerated) {
    const auto kmer_index_entry = "5 9 @7 19";
    const std::vector<std::string> &parts = split(kmer_index_entry, "|");

    const auto &result = parse_variant_site_path(parts[0]);
    VariantSitePath expected = {
            VariantSite {5, 9},
            VariantSite {7, 19},
    };
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenSitesTrailingAt_TrailingAtIgnored) {
    const auto kmer_index_entry = "5 9 @7 19 @";
    const std::vector<std::string> &parts = split(kmer_index_entry, "|");

    const auto &result = parse_variant_site_path(parts[0]);
    VariantSitePath expected = {
            VariantSite {5, 9},
            VariantSite {7, 19},
    };
    EXPECT_EQ(result, expected);
}


class IndexKmers : public ::testing::Test {

protected:
    std::string prg_fpath;

    virtual void SetUp() {
        boost::uuids::uuid uuid = boost::uuids::random_generator()();
        const auto uuid_str = boost::lexical_cast<std::string>(uuid);
        prg_fpath = "./prg_" + uuid_str;
    }

    virtual void TearDown() {
        std::remove(prg_fpath.c_str());
    }

    FM_Index fm_index_from_raw_prg(const std::string &prg_raw) {
        std::vector<uint64_t> prg = encode_prg(prg_raw);
        dump_encoded_prg(prg, prg_fpath);
        FM_Index fm_index;
        // TODO: constructing from memory with sdsl::construct_im appends 0 which corrupts
        sdsl::construct(fm_index, prg_fpath, 8);
        return fm_index;
    }

    PRG_Info generate_prg_info(const std::string &prg_raw) {
        PRG_Info prg_info;
        prg_info.fm_index = fm_index_from_raw_prg(prg_raw);
        prg_info.dna_rank = calculate_ranks(prg_info.fm_index);
        prg_info.allele_mask = generate_allele_mask(prg_raw);
        prg_info.max_alphabet_num = max_alphabet_num(prg_raw);
        return prg_info;
    }

};


TEST_F(IndexKmers, KmerCrossesVariantRegion_KmerNotInNonVariantRegionSet) {
    const std::string prg_raw = "aca5g6t5gcatt";
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("atgca");
    const int kmer_size = 5;
    Patterns kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.non_site_crossing_kmers;
    NonSiteCrossingKmers expected = {};
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, KmerInNonVariantRegion_KmerIncludedInNonVarKmerSet) {
    const std::string prg_raw = "aca5g6t5gcatt";
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("gcatt");
    const int kmer_size = 5;
    Patterns kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.non_site_crossing_kmers;
    NonSiteCrossingKmers expected = {kmer};
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


TEST_F(IndexKmers, KmerCrossesSecondAllele_VariantRegionRecordedInSites) {
    const std::string prg_raw = "aca5g6t5gctc";
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("atgct");
    const int kmer_size = 5;
    Patterns kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.variant_site_paths_map[kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 2}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, KmerCrossesFirstAllele_VariantRegionRecordedInSites) {
    const std::string prg_raw = "aca5g6t5gcatt";
    const auto prg_info = generate_prg_info(prg_raw);

    auto kmer = encode_dna_bases("aggca");
    const int kmer_size = 5;
    Patterns kmers = {kmer};

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.variant_site_paths_map[kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 1}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, BothKmersOverlapVariantSiteAlleles_CorrectSearchResults) {
    const std::string prg_raw = "aca5g6c5tatt";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto kmer_suffix_diff = encode_dna_bases("ac");
    auto second_full_kmer = encode_dna_bases("actat");
    Patterns kmers = {
            first_full_kmer,
            kmer_suffix_diff
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.variant_site_paths_map;
    KmerVariantSitePaths expected = {
            {first_full_kmer,
                    VariantSitePaths {
                            VariantSitePath {
                                    VariantSite {5, 1}
                            }
                    }
            },

            {second_full_kmer,
                    VariantSitePaths {
                            VariantSitePath {
                                    VariantSite {5, 2}
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, OneKmersOverlapsVariantSiteAllele_CorrectSearchResults) {
    const std::string prg_raw = "aca5g6c5tatt";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 5;
    auto first_full_kmer = encode_dna_bases("agtat");
    auto kmer_suffix_diff = encode_dna_bases("aa");
    auto second_full_kmer = encode_dna_bases("aatat");
    Patterns kmers = {
            first_full_kmer,
            kmer_suffix_diff
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.variant_site_paths_map[first_full_kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 1}
            }
    };
    EXPECT_EQ(result, expected);

    result = kmer_index.variant_site_paths_map[second_full_kmer];
    expected = {};
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, ThreeKmersOverlapSiteThreeAllele_CorrectSearchResults) {
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

    auto &result = kmer_index.variant_site_paths_map[first_full_kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 1}
            }
    };
    EXPECT_EQ(result, expected);

    result = kmer_index.variant_site_paths_map[second_full_kmer];
    expected = {
            VariantSitePath {
                    VariantSite {5, 2}
            }
    };
    EXPECT_EQ(result, expected);

    result = kmer_index.variant_site_paths_map[third_full_kmer];
    expected = {
            VariantSitePath {
                    VariantSite {5, 3}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, ThreeKmersOneMissMatch_CorrectSearchResults) {
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

    auto &result = kmer_index.variant_site_paths_map[first_full_kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 1}
            }
    };
    EXPECT_EQ(result, expected);

    result = kmer_index.variant_site_paths_map[second_full_kmer];
    expected = {
            VariantSitePath {
                    VariantSite {5, 2}
            }
    };
    EXPECT_EQ(result, expected);

    result = kmer_index.variant_site_paths_map[third_full_kmer];
    expected = {};
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, OneKmerStartsAtAllele_SiteFound) {
    const std::string prg_raw = "aca5g6c6a5tatt";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("gtat");
    Patterns kmers = {
            encode_dna_bases("gtat"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.variant_site_paths_map[first_full_kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 1}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, TwoKmersStartAtAllele_SitesFound) {
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

    auto &result = kmer_index.variant_site_paths_map[first_full_kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 1}
            }
    };
    EXPECT_EQ(result, expected);

    result = kmer_index.variant_site_paths_map[second_full_kmer];
    expected = {
            VariantSitePath {
                    VariantSite {5, 2}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, KmerEndingInAllele_SingleSiteFound) {
    const std::string prg_raw = "aca5g6c5t";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("acag");
    Patterns kmers = {
            encode_dna_bases("acag"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.variant_site_paths_map[first_full_kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 1}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, TwoKmersEndingInAlleles_TwoSingleSitesFound) {
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

    auto &result = kmer_index.variant_site_paths_map[first_full_kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 1}
            }
    };
    EXPECT_EQ(result, expected);

    result = kmer_index.variant_site_paths_map[second_full_kmer];
    expected = {
            VariantSitePath {
                    VariantSite {5, 2}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST_F(IndexKmers, KmerStartingInSiteAndEndInAnotherSite_CorrectVariantSitePath) {
    const std::string prg_raw = "aca5g6c5tt7a8c7gg";
    const auto prg_info = generate_prg_info(prg_raw);

    const int kmer_size = 4;
    auto first_full_kmer = encode_dna_bases("ctta");
    Patterns kmers = {
            encode_dna_bases("ctta"),
    };

    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

    auto &result = kmer_index.variant_site_paths_map[first_full_kmer];
    VariantSitePaths expected = {
            VariantSitePath {
                    VariantSite {5, 2},
                    VariantSite {7, 1}
            }
    };
    EXPECT_EQ(result, expected);
}
