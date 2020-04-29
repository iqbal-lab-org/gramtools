#include <filesystem>

#include "gtest/gtest.h"

#include "prg/linearised_prg.hpp"

/************************/
/* Conversion utilities */
/************************/

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

/********************/
/* PRG_String class */
/********************/

namespace fs = std::filesystem;
auto const test_data_dir = fs::path(__FILE__).parent_path().parent_path() / "test_data";

TEST(PRGString, Load_from_File){
    /*
     * The tested file is the binary output of running `make_prg` on the following MSA:
                             ">R1\n"
                             "AAAAAAAAA\n"
                             ">R2\n"
                             "AATAAAAAA\n"
                             ">R3\n"
                             "AAAAATAAA\n"
                             ">R4\n"
                             "TTTTTTTTT\n"
                             ">R5\n"
                             "TTATTTTTT\n"
                             ">R6\n"
                             "TTTTTATTT\n";
     */
    fs::path path(test_data_dir / "twoSegregatingClasses.fasta.max_nest10.min_match1.bin");
    PRG_String l = PRG_String(path.generic_string());
    std::string expected{"[AA[A,T]AA[A,T]AAA,TT[A,T]TT[A,T]TTT]"};
    auto res = ints_to_prg_string(l.get_PRG_string());
    EXPECT_EQ(expected, res);
}

class PRGString_WriteAndRead : public ::testing::Test {
protected:
    void SetUp() {
        fname = "@pstring_out";
        std::string prg_string{"A[A,C]T[GGG,G]C"};
        expected_markers = prg_string_to_ints(prg_string);
        p = PRG_String(expected_markers);
    }

    void TearDown() {
        if (remove(fname.c_str()) != 0){
            std::cerr << "Could not delete the built file " << fname;
            exit(1);
        }
    }
    std::string fname;
    marker_vec expected_markers;
    PRG_String p;
};

TEST_F(PRGString_WriteAndRead, WriteAndReadLittleEndian){
    // Little Endian should be the default for read and write
    p.write(fname);

    // Load it into another object
    PRG_String p2{fname};
    EXPECT_EQ(p2.get_endianness(), endianness::little);
    EXPECT_EQ(expected_markers, p2.get_PRG_string());
}

TEST_F(PRGString_WriteAndRead, WriteAndReadBigEndian){
    p.write(fname, endianness::big);

    // Load it into another object
    PRG_String p2{fname, endianness::big};
    EXPECT_EQ(p2.get_endianness(), endianness::big);
    EXPECT_EQ(expected_markers, p2.get_PRG_string());
}

TEST(PRGString, ExitPoint_MapPositions){
    marker_vec t{5,1,6,2,7,1,8,3,8,6}; // Ie: "[A,C[A,T]]"
    PRG_String l = PRG_String(t);
    std::unordered_map<Marker,int> expected_end_positions{
            {6 ,9},
            {8, 8}
    };
    EXPECT_EQ(expected_end_positions, l.get_end_positions());
}
