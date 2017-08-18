#include "gtest/gtest.h"
#include "kmers.hpp"


TEST(GeneratePrecalc, GivenDataForSinglePrecalcEntry_CorrectDumpRowGenerated) {
    const Kmer kmer = {1, 2, 3, 4};
    const KmersRef kmers_in_ref = {kmer};

    Site first_site = {
            VariantSite(5, {9, 8, 7}),
            VariantSite(7, {19, 18, 17})
    };
    Site second_site = {
            VariantSite(9, {29, 28, 27}),
            VariantSite(11, {39, 38, 37})
    };
    Sites sites = {first_site, second_site};
    const KmerSites kmer_sites = {
            {kmer, sites},
    };

    const SA_Intervals sa_intervals = {
            {123, 456},
            {789, 424}
    };

    const auto result = dump_kmer_precalc_entry(kmer,
                                                sa_intervals,
                                                kmers_in_ref,
                                                kmer_sites);
    const auto expected = "1 2 3 4|1|123 456 789 424||5 9 8 7 @7 19 18 17 @|9 29 28 27 @11 39 38 37 @|";
    EXPECT_EQ(result, expected);
}


TEST(GeneratePrecalc, GivenSites_DumpSitesCorrectly) {
    Site first_site = {
            VariantSite(5, {9, 8, 7}),
            VariantSite(7, {19, 18, 17})
    };
    Site second_site = {
            VariantSite(9, {29, 28, 27}),
            VariantSite(11, {39, 38, 37})
    };
    Sites sites = {first_site, second_site};

    const Kmer kmer = {1, 2, 3, 4};
    const KmerSites kmer_sites = {
            {kmer, sites},
    };
    const auto result = dump_sites(kmer, kmer_sites);
    const auto expected = "5 9 8 7 @7 19 18 17 @|9 29 28 27 @11 39 38 37 @|";
    EXPECT_EQ(result, expected);
}


TEST(GeneratePrecalc, GivenSaIntervals_DumpSaIntervalsStringCorrectly) {
    SA_Intervals sa_intervals = {
            SA_Interval(1, 2),
            SA_Interval(3, 4)
    };
    const auto result = dump_sa_intervals(sa_intervals);
    const auto expected = "1 2 3 4";
    EXPECT_EQ(result, expected);
}


TEST(GeneratePrecalc, GivenKmer_DumpKmerStringCorrectly) {
    const Kmer kmer = {1, 2, 3, 4};
    const auto result = dump_kmer(kmer);
    const auto expected = "1 2 3 4";
    EXPECT_EQ(result, expected);
}


TEST(GeneratePrecalc, GivenDnaString_DnaBasesEncodedCorrectly) {
    const auto dna_str = "AAACCCGGGTTTACGT";
    const auto result = encode_dna_bases(dna_str);
    const std::vector<uint8_t> expected = {
            1, 1, 1,
            2, 2, 2,
            3, 3, 3,
            4, 4, 4,
            1, 2, 3, 4,
    };
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenEncodedKmerString_CorrectlyParsed) {
    const auto encoded_kmer_str = "3 4 2 1 1 3 1 1 2";
    const auto result = parse_encoded_kmer(encoded_kmer_str);
    const Kmer expected = {3, 4, 2, 1, 1, 3, 1, 1, 2};
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenSaIntervalsString_CorrectlyParsed) {
    const auto full_sa_intervals_str = "352511 352512 352648 352649 352648 352649";
    const auto result = parse_sa_intervals(full_sa_intervals_str);

    SA_Intervals expected {
            {352511, 352512},
            {352648, 352649},
            {352648, 352649},
    };
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenTwoSites_CorrectSiteStructGenerated) {
    Site expected = {
            VariantSite(5, {9, 8, 7}),
            VariantSite(7, {19, 18, 17})
    };

    const auto precalc_kmer_entry = "5 9 8 7 @7 19 18 17";
    const std::vector<std::string> &parts = split(precalc_kmer_entry, "|");
    const auto &result = parse_site(parts[0]);
    EXPECT_EQ(result, expected);
}


TEST(ParsePrecalc, GivenSitesTrailingAt_TrailingAtIgnored) {
    Site expected = {
            VariantSite(5, {9, 8, 7}),
            VariantSite(7, {19, 18, 17})
    };

    const auto precalc_kmer_entry = "5 9 8 7 @7 19 18 17 @";
    const std::vector<std::string> &parts = split(precalc_kmer_entry, "|");
    const auto &result = parse_site(parts[0]);
    EXPECT_EQ(result, expected);
}
