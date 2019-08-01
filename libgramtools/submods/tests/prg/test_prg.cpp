#include "gtest/gtest.h"
#include "src_common/generate_prg.hpp"

#include "prg/masks.hpp"


using namespace gram;


TEST(GetMaxAlphabetNum, GivenPrg_CorrectMaxAlphabetNum) {
    auto prg_raw = "a5g6t6cccc11g12tttt12";
    auto prg_info = generate_prg_info(prg_raw);
    auto result = get_max_alphabet_num(prg_info.encoded_prg);
    auto expected = 12;
    EXPECT_EQ(result, expected);
}


TEST(GetMaxAlphabetNum, PrgWithVariantSite_LargestSiteMarkerAsMaxAlphabet) {
    auto prg_raw = "a13g14t14tt";
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
    auto prg_raw = "a5g6t6cc11g12tt12";
    auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_sites_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {0, 0, 5, 0, 5, 0, 0, 0, 0, 11, 0, 11, 11, 0};
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(GenerateSitesMask, SingleVariantSiteTwoAlleles_CorrectSitesMask) {
    const std::string prg_raw = "a5g6t6c";
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_sites_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0, 0, 5, 0, 5, 0, 0
    };
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(GenerateSitesMask, TwoVariantSites_CorrectSitesMask) {
    const std::string prg_raw = "a5g6t6cc7g8tt8aa8";
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
