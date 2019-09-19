#include "gtest/gtest.h"
#include "common/utils.hpp"

using namespace gram;

TEST(ReverseComplementRead, GivenRead_ReverseComplementReadReturned) {
    gram::Pattern read = {1, 2, 1, 3, 4};
    auto result = gram::reverse_complement_read(read);
    gram::Pattern expected = {1, 2, 4, 3, 4};
    EXPECT_EQ(result, expected);
}

TEST(PRG_Conversion, string_to_ints1){
    std::string prg_string = "[A,C[A,T]]";
    std::vector<Marker> expected = {5,1,6,2,7,1,8,4,8,6};
    auto res = prg_string_to_ints(prg_string);
    EXPECT_EQ(res, expected);
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