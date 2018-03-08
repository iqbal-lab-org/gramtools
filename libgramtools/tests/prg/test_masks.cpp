#include "gtest/gtest.h"

#include "../test_utils.hpp"
#include "prg/masks.hpp"


TEST(LoadAlleleMask, GivenComplexAlleleMask_SaveAndLoadFromFileCorrectly) {
    auto prg_raw = "a5g6ttt5cc7aa8t7a";
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
    auto prg_raw = "a5g6ttt5cc7aa8t7a";
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
    const std::string prg_raw = "a5g6t5c";
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
    const std::string prg_raw = "a5g6t6aa5c";
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
    const std::string prg_raw = "a5g6t5cc7aa8g7a";
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
    const std::string prg_raw = "a13g14t13tt";
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