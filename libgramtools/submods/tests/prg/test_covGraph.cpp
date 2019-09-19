#include "gtest/gtest.h"
#include "build/coverage_graph.hpp"

using namespace gram;


/*
 * -----------------------
 * `PRG String` tests
 * -----------------------
 */

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

/*
 * -----------------------
 * `cov_Graph_Builder` tests
 * NOTE: the best way to understand these tests is to draw the DAG corresponding to the PRG String being tested,
 * labeling nodes with their expected attributes (eg position, site/allele ID).
 *
 * Use test fixtures: single data, multiple tests
 * -----------------------
 */
class cov_G_Builder_nested : public ::testing::Test {
protected:
    void SetUp() {
        std::string prg_string{"[A,AA,A[A,C]A]C[AC,C]G"}; // A simple nested string
        //idx:                    0    5     11     18
        marker_vec v = prg_string_to_ints(prg_string);
        PRG_String p{v};
        p.process();
        c = cov_Graph_Builder{p};
    }
    cov_Graph_Builder c;
};

TEST_F(cov_G_Builder_nested, FindMarkerTypes){
    //"[A,AA,A[A,C]A]C[AC,C]G"
    using mt = marker_type;
   marker_type res;
   std::array<int, 5> positions{0, 2, 4, 11, 13};
   std::array<marker_type, 5> expected{mt::site_entry, mt::allele_end, mt::sequence, mt::site_end, mt::site_end};

   for (int i = 0; i < positions.size(); ++i){
       res = c.find_marker_type(positions[i]);
       EXPECT_EQ(res,expected[i]);
   }
}

TEST_F(cov_G_Builder_nested, Bubbles_Positions){
    //"[A,AA,A[A,C]A]C[AC,C]G"
    c.run();
    // The positions are not INDICES in the PRG string; they are the positions in the multiple-sequence alignment
    // giving rise to it. Draw the graph of the PRG string and take the LONGEST allele positions to obtain them.
    std::vector start_positions{0, 1, 4}; // Take the positions of the '[' in left to right order in PRG string
    std::vector end_positions{3, 2, 6}; // To know these, take the ']' position corresponding to the '[' start position.

    // The bubbles are in 'topological order': the highest start position is processed FIRST
    // So we test the vector indices high to low
    int bubble_number = start_positions.size() - 1;
    for (auto& bubble : c.bubble_map){
        EXPECT_EQ(bubble.first->get_pos(), start_positions[bubble_number]);
        EXPECT_EQ(bubble.second->get_pos(), end_positions[bubble_number]);
        bubble_number--;
    }
}

TEST_F(cov_G_Builder_nested, ParentalMap){
    //"[A,AA,A[A,C]A]C[AC,C]G"
    c.run();
    // Expecting to find a single entry, for the single nested site, pointing to siteID 5 & alleleID 3.
    parental_map expected {
            {7 , VariantLocus{5, 3} }
    };
    EXPECT_EQ(c.par_map, expected);
}
