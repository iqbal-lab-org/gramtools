#include "gtest/gtest.h"
#include "utils.hpp"
#include "prg.hpp"


TEST(GenerateAlleleMask, SingleVariantSite_CorrectAlleleMask) {
    const std::string prg_raw = "a5g6t5c";
    auto result = generate_allele_mask(prg_raw);
    std::vector<AlleleId> expected = {
            0,
            0, 1,
            0, 2, 0,
            0
    };
    EXPECT_EQ(result, expected);
}


TEST(GenerateAlleleMask, SingleVariantSiteThreeAlleles_CorrectAlleleMask) {
    const std::string prg_raw = "a5g6t6aa5c";
    auto result = generate_allele_mask(prg_raw);
    std::vector<AlleleId> expected = {
            0,
            0, 1,
            0, 2,
            0, 3, 3, 0,
            0
    };
    EXPECT_EQ(result, expected);
}


TEST(GenerateAlleleMask, TwoVariantSites_CorrectAlleleMask) {
    const std::string prg_raw = "a5g6t5cc7aa8g7a";
    auto result = generate_allele_mask(prg_raw);
    std::vector<AlleleId> expected = {
            0,
            0, 1,
            0, 2, 0,
            0, 0,
            0, 1, 1,
            0, 2, 0,
            0,
    };
    EXPECT_EQ(result, expected);
}


TEST(GenerateAlleleMask, DoubleDigitMarker_CorrectAlleleMask) {
    const std::string prg_raw = "a13g14t13tt";
    auto result = generate_allele_mask(prg_raw);
    std::vector<AlleleId> expected = {
            0,
            0, 1,
            0, 2, 0,
            0, 0,
    };
    EXPECT_EQ(result, expected);
}

TEST(MaxAlphabetNum, PrgWithVariantSite_LargestSiteMarkerAsMaxAlphabet) {
    const std::string prg_raw = "a13g14t13tt";
    auto result = max_alphabet_num(prg_raw);
    uint64_t expected = 14;
    EXPECT_EQ(result, expected);
}

TEST(MaxAlphabetNum, SingleCharPrg_CorrectBaseEncodingAsMaxAlphabet) {
    const std::string prg_raw = "c";
    auto result = max_alphabet_num(prg_raw);
    uint64_t expected = 2;
    EXPECT_EQ(result, expected);
}
