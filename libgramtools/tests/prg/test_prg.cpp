#include "gtest/gtest.h"
#include "../test_utils.hpp"

#include "prg/masks.hpp"


using namespace gram;


TEST(GetMaxAlphabetNum, GivenPrg_CorrectMaxAlphabetNum) {
    auto prg_raw = "a5g6t5cccc11g12tttt11";
    auto prg_info = generate_prg_info(prg_raw);
    auto result = get_max_alphabet_num(prg_info.encoded_prg);
    auto expected = 12;
    EXPECT_EQ(result, expected);
}


TEST(GetMaxAlphabetNum, PrgWithVariantSite_LargestSiteMarkerAsMaxAlphabet) {
    auto prg_raw = "a13g14t13tt";
    auto prg_info = generate_prg_info(prg_raw);
    auto result = get_max_alphabet_num(prg_info.encoded_prg);
    uint64_t expected = 14;
    EXPECT_EQ(result, expected);
}


TEST(GetMaxAlphabetNum, SingleCharPrg_CorrectBaseEncodingAsMaxAlphabet) {
    auto prg_raw = "c";
    auto prg_info = generate_prg_info(prg_raw);
    auto result = get_max_alphabet_num(prg_info.encoded_prg);
    uint64_t expected = 2;
    EXPECT_EQ(result, expected);
}


TEST(GenerateSitesMask, GivenMultiSitePrg_CorrectSitesMask) {
    auto prg_raw = "a5g6t5cc11g12tt11";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_sites_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {0, 0, 5, 0, 5, 0, 0, 0, 0, 11, 0, 11, 11, 0};
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(GenerateSitesMask, SingleVariantSiteTwoAlleles_CorrectSitesMask) {
    const std::string prg_raw = "a5g6t5c";
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_sites_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0, 0, 5, 0, 5, 0, 0
    };
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(GenerateSitesMask, TwoVariantSites_CorrectSitesMask) {
    const std::string prg_raw = "a5g6t5cc7g8tt8aa7";
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_sites_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0,
            0, 5, 0, 5, 0,
            0, 0,
            0, 7, 0, 7, 7, 0, 7, 7, 0
    };
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}
