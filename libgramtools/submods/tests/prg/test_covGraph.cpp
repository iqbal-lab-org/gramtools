#include "gtest/gtest.h"
#include "build/load_PRG_string.hpp"

TEST(Load_PRGString, l){
    std::string fin = "/home/brice/Desktop/git_repos/prg_string_construction/python_make_prg/test_data/twoSegregatingClasses.fasta.max_nest10.min_match7.bin";
    PRG_String l = PRG_String(fin);
    for (auto& s : l.get_PRG_string()) std::cout << s;
}

TEST(Load_PRGString, fromIntVector){
   int_vec t{5,1,6,2,5};
    PRG_String l = PRG_String(t);
    //EXPECT_EQ(t, nullptr);
    int_vec expected{5,1,6,2,5};
    EXPECT_EQ(expected, l.get_PRG_string());

    l.process();
    expected = {5,1,6,2,6};
    EXPECT_EQ(expected, l.get_PRG_string());
}
