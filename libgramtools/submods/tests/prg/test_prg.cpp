#include "gtest/gtest.h"
#include "src_common/generate_prg.hpp"

#include "kmer_index/masks.hpp"


using namespace gram;

TEST(GetNumVarSites, NoSites) {
    auto prg_raw = encode_prg("c");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = prg_info.num_variant_sites;
    uint64_t expected = 0;
    EXPECT_EQ(result, expected);
}

TEST(GetNumVarSites, UnNestedPrgString) {
    auto prg_raw = encode_prg("a5g6t6cccc11g12tttt12");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = prg_info.num_variant_sites;
    auto expected = 2;
    EXPECT_EQ(result, expected);
}


TEST(GetNumVarSites, Nested_PrgString) {
    auto prg_raw = prg_string_to_ints("[[A,C,G]A,T]T[,C][GA,CT]");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = prg_info.num_variant_sites;
    uint64_t expected = 4;
    EXPECT_EQ(result, expected);
}

TEST(GenerateSitesMask, GivenMultiSitePrg_CorrectSitesMask) {
    auto prg_raw = encode_prg("a5g6t6cc11g12tt12");
    auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_sites_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {0, 0, 5, 0, 5, 0, 0, 0, 0, 11, 0, 11, 11, 0};
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(GenerateSitesMask, SingleVariantSiteTwoAlleles_CorrectSitesMask) {
    const auto prg_raw = encode_prg("a5g6t6c");
    auto prg_info = generate_prg_info(prg_raw);
    auto result = generate_sites_mask(prg_info.encoded_prg);
    sdsl::int_vector<> expected = {
            0, 0, 5, 0, 5, 0, 0
    };
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(GenerateSitesMask, TwoVariantSites_CorrectSitesMask) {
    const auto prg_raw = encode_prg("a5g6t6cc7g8tt8aa8");
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
