#include "gtest/gtest.h"
#include "build/load_PRG_string.hpp"

using namespace gram;

TEST(PRGString, Load_from_File){
    std::string fin = "/home/brice/Desktop/git_repos/prg_string_construction/python_make_prg/test_data/twoSegregatingClasses.fasta.max_nest10.min_match7.bin";
    PRG_String l = PRG_String(fin);
    for (auto& s : l.get_PRG_string()) std::cout << s;
}

TEST(PRGString, Load_fromIntVector){
   marker_vec t{5,1,6,2,5};
    PRG_String l = PRG_String(t);
    marker_vec expected{5,1,6,2,5};
    EXPECT_EQ(expected, l.get_PRG_string());
}

TEST(PRGString, ExitPoint_ConvertOddToEven){
    marker_vec t{5,1,6,2,5};
    PRG_String l = PRG_String(t);
    l.process();
    EXPECT_EQ(true, l.odd_site_end_found);
    // The vector should now have even site marker exit points
    marker_vec expected = {5,1,6,2,6};
    EXPECT_EQ(expected, l.get_PRG_string());
}

TEST(PRGString, ExitPoint_MapPositions){
    marker_vec t{5,1,6,2,7,1,8,3,8,6}; // Ie: "[A,C[A,T]]"
    PRG_String l = PRG_String(t);
    l.process();
    std::unordered_map<Marker,int> expected_end_positions{
            {6 ,9},
            {8, 8}
    };
    EXPECT_EQ(expected_end_positions, l.get_end_positions());
}
