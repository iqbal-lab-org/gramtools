#include "gtest/gtest.h"
#include "submod_resources.hpp"

#include "build/kmer_index/masks.hpp"


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


TEST(PRG_Conversion, string_to_ints1){
    std::string prg_string = "[A,C[A,T]]";
    std::vector<Marker> expected = {5,1,6,2,7,1,8,4,8,6};
    auto res = prg_string_to_ints(prg_string);
    EXPECT_EQ(res, expected);
}


TEST(PRG_Conversion, StringWithInvalidCharPassed_ProgramExits){
    std::string prg_string = "5A5";
    ASSERT_DEATH(prg_string_to_ints(prg_string),  ".*not a nucleotide*.");
}

TEST(PRG_Conversion, ints_to_string){
    std::vector<Marker> int_vec = {5,1,6,2,7,1,8,4,8,6};
    std::string expected = "[A,C[A,T]]";
    auto res = ints_to_prg_string(int_vec);
    EXPECT_EQ(res, expected);
}

TEST(PRG_Conversion, string_to_ints2){
    std::string prg_string = "[AAA,,A[CCC,CC,C]]G";
    std::vector<Marker> expected = {5,1,1,1,6,6,1,7,2,2,2,8,2,2,8,2,8,6,3};
    auto res = prg_string_to_ints(prg_string);
    EXPECT_EQ(res, expected);
}

TEST(PRG_Conversion, string_to_ints3) {
    std::string prg_string = "[A,AA,A[A,C]A]C[A,C]";
    std::vector<Marker> expected = {5,1,6,1,1,6,1,7,1,8,2,8,1,6,2,9,1,10,2,10};
    auto res = prg_string_to_ints(prg_string);
    EXPECT_EQ(res, expected);
}

/**
 * Here I want to highlight that the initial site numbering gets lost by int to string conversion
 * if the initial site numbering does not obey: 'sites entered first have smaller site IDs'
 */
TEST(PRG_Conversion, ints_to_string_to_ints){
    std::vector<Marker> int_vec = {7,1,8,2,5,1,6,4,6,8};
    std::string expected_string = "[A,C[A,T]]";
    auto res1 = ints_to_prg_string(int_vec);
    EXPECT_EQ(res1, expected_string);

    std::vector<Marker> expected_vec = {5,1,6,2,7,1,8,4,8,6};
    auto res2 = prg_string_to_ints(expected_string);
    EXPECT_EQ(res2, expected_vec);
}

TEST(BuildChildMap, GivenParentalMap_CorrectChildMap){
    // Site 5 has two sites nested in haplogroup 1, and one in haplogroup 2.
    // Note: parental_map / quasimap stores allele haplogroups as 1-based,
    // but child_map moves them to 0-based (consistent with infer).
    parental_map par_map{
            {7, VariantLocus{5, FIRST_ALLELE}},
            {9, VariantLocus{5, FIRST_ALLELE}},
            {11, VariantLocus{5, FIRST_ALLELE + 1}},
            {15, VariantLocus{13, FIRST_ALLELE + 2}},
    };

    auto result = build_child_map(par_map);
    // Sort the internal vectors to be independent from parental hash map element ordering
    for (auto& entry : result){
        for (auto& entry2 : entry.second){
            std::sort(entry2.second.begin(), entry2.second.end());
        }
    }
    child_map expected{
            {
                    5, haplo_map{
                    {0, marker_vec{7, 9} },
                    {1, marker_vec{11} },
            }
            },
            {
                    13, haplo_map{
                    {2, marker_vec{15}}
            }
            }
    };

    EXPECT_EQ(result, expected);
}
