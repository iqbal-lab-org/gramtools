#include "gtest/gtest.h"
#include "../test_utils.hpp"


TEST(GetMaxAlphabetNum, GivenPrg_CorrectMaxAlphabetNum) {
    auto prg_raw = "a5g6t5cccc11g12tttt11";
    auto prg_info = generate_prg_info(prg_raw);
    auto result = get_max_alphabet_num(prg_info.encoded_prg);
    auto expected = 12;
    EXPECT_EQ(result, expected);
}
