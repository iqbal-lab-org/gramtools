#include "gtest/gtest.h"

#include "process_prg.hpp"
#include "map.hpp"


TEST(ConstructVariantSitesMask, CorrectBitsSet) {
    std::string test_fpath = "./test_cases/one_snp.txt";
    FM_Index fm_index = construct_fm_index(test_fpath,
                                           "int_alphabet_file",
                                           "memory_log_file",
                                           "csa_file", true);

    sdsl::bit_vector variant_mask = construct_variant_sites_mask(fm_index);
    sdsl::bit_vector expected = {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0};
    EXPECT_EQ(variant_mask, expected);
}