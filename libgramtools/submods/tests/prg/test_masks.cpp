#include "gtest/gtest.h"

#include "src_common/generate_prg.hpp"
#include "kmer_index/masks.hpp"


using namespace gram;


TEST(LoadAlleleMask, GivenComplexAlleleMask_SaveAndLoadFromFileCorrectly) {
    auto prg_raw = encode_prg("a5g6ttt6cc7aa8t8a");
    auto prg_info = generate_prg_info(prg_raw);
    auto allele_mask = generate_allele_mask(prg_info.encoded_prg);

    Parameters parameters = {};
    parameters.allele_mask_fpath = "@allele_mask";
    sdsl::store_to_file(allele_mask, parameters.allele_mask_fpath);

    auto result = load_allele_mask(parameters);
    sdsl::int_vector<> expected = {
            0,
            0, 1, 0, 2, 2, 2, 0,
            0, 0,
            0, 1, 1, 0, 2, 0,
            0
    };
    for (auto i = 0; i < result.size(); ++i)
        EXPECT_EQ(result[i], expected[i]);
}


TEST(GenerateAlleleMask, GivenMultipleSitesAndAlleles_CorrectAlleleMask) {
    auto prg_raw = encode_prg("a5g6ttt6cc7aa8t8a");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_allele_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0,
            0, 1, 0, 2, 2, 2, 0,
            0, 0,
            0, 1, 1, 0, 2, 0,
            0
    };
    for (auto i = 0; i < result.size(); ++i)
        EXPECT_EQ(result[i], expected[i]);
}


TEST(GenerateAlleleMask, SingleVariantSite_CorrectAlleleMask) {
    const auto prg_raw = encode_prg("a5g6t6c");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_allele_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0,
            0, 1,
            0, 2, 0,
            0
    };
    for (auto i = 0; i < result.size(); ++i)
        EXPECT_EQ(result[i], expected[i]);
}


TEST(GenerateAlleleMask, SingleVariantSiteThreeAlleles_CorrectAlleleMask) {
    const auto prg_raw = encode_prg("a5g6t6aa6c");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_allele_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0,
            0, 1,
            0, 2,
            0, 3, 3, 0,
            0
    };
    for (auto i = 0; i < result.size(); ++i)
        EXPECT_EQ(result[i], expected[i]);
}


TEST(GenerateAlleleMask, TwoVariantSites_CorrectAlleleMask) {
    const auto prg_raw = encode_prg("a5g6t6cc7aa8g8a");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_allele_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0,
            0, 1,
            0, 2, 0,
            0, 0,
            0, 1, 1,
            0, 2, 0,
            0,
    };
    for (auto i = 0; i < result.size(); ++i)
        EXPECT_EQ(result[i], expected[i]);
}


TEST(GenerateAlleleMask, DoubleDigitMarker_CorrectAlleleMask) {
    const auto prg_raw = encode_prg("a13g14t14tt");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_allele_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0,
            0, 1,
            0, 2, 0,
            0, 0,
    };
    for (auto i = 0; i < result.size(); ++i)
        EXPECT_EQ(result[i], expected[i]);
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

TEST(GenerateBWTMask, rankQueries){
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
