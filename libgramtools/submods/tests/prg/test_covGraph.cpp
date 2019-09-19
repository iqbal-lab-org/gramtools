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
        marker_vec v = prg_string_to_ints(prg_string);
        PRG_String p{v};
        p.process();
        c = cov_Graph_Builder{p};
    }
    cov_Graph_Builder c;
};

// Test that marker typing is correct
TEST_F(cov_G_Builder_nested, FindMarkerTypes){
    //                       "[A,AA,A[A,C]A]C[AC,C]G"
    //idx:                    0    5     11     18
    using mt = marker_type;
   marker_type res;
   std::array<int, 5> positions{0, 2, 4, 11, 13};
   std::array<marker_type, 5> expected{mt::site_entry, mt::allele_end, mt::sequence, mt::site_end, mt::site_end};

   for (int i = 0; i < positions.size(); ++i){
       res = c.find_marker_type(positions[i]);
       EXPECT_EQ(res,expected[i]);
   }
}


// Test that the parental map is correct
TEST_F(cov_G_Builder_nested, ParentalMap){
    //"[A,AA,A[A,C]A]C[AC,C]G"
    c.run();
    // Expecting to find a single entry, for the single nested site, pointing to siteID 5 & alleleID 3.
    parental_map expected {
            {7 , VariantLocus{5, 3} }
    };
    EXPECT_EQ(c.par_map, expected);
}

// Test that the node site & allele IDs are correct
TEST_F(cov_G_Builder_nested, SiteAndAllele_IDs){
    //"[A,AA,A[A,C]A]C[AC,C]G"
    c.run();
    auto const& rand_access = c.random_access;
    std::vector<VariantLocus> expected{
            {5, 0}, {5, 1}, {5, 0}, {5, 2}, {5, 2},
            {5, 0}, {5, 3}, {7, 0}, {7, 1}, {7, 0},
            {7, 2}, {7, 0}, {5, 3}, {5, 0}, {0, 0},
            {9, 0}, {9, 1}, {9, 1}, {9, 0}, {9,2},
            {9, 0}, {0, 0}
    };
    std::vector<VariantLocus> res{rand_access.size(), VariantLocus()};
    int pos = 0;
    for (auto const&s : rand_access){
       res[pos].first = s.node->get_site();
        res[pos].second = s.node->get_allele();
        pos++;
    }

    EXPECT_EQ(res, expected);
}

// Test that the node positions are correct
TEST_F(cov_G_Builder_nested, SequencePositions) {
    //"[A,AA,A[A,C]A]C[AC,C]G"
    c.run();
    auto const &rand_access = c.random_access;
    // The positions are not INDICES in the PRG string; they are the positions in the multiple-sequence alignment
    // giving rise to it. Draw the graph of the PRG string and take the LONGEST allele positions to obtain them.
    std::vector<seqPos> expected{
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, // First site exit point here
                2, 3, 3, 4, 4, 4, 4, 4, 6, 6 };
    std::vector<seqPos> res(expected.size(), 0);
    int pos = 0;
    for (auto const &s : rand_access) {
        res[pos++] = s.node->get_pos();
    }
    EXPECT_EQ(res, expected);
}

// Test that bubble entry and exit points are correctly identified
TEST_F(cov_G_Builder_nested, Bubbles){
    //"[A,AA,A[A,C]A]C[AC,C]G
    using i_v = std::vector<int>; // Stores indexes into PRG string
    c.run();
    auto const &rand_access = c.random_access;
    // Note: allele separators (",") point to the site entry node, so we expect them here
    i_v expected_site_entry_points{
       0, 2, 5, 7, 9, 15, 18 };
    i_v expected_site_exit_points{
            11, 13, 20 };
    i_v res_entries, res_exits;

    int pos = -1;
    for (auto const &s : rand_access) {
        pos++;
        Marker site_ID = s.node->get_site();
        try{
            bool is_site_entry = c.bubble_starts.at(site_ID) == s.node;
            bool is_site_exit = c.bubble_ends.at(site_ID) == s.node;
            EXPECT_FALSE(is_site_entry & is_site_exit); // They should not both be true
            if (is_site_entry){
                EXPECT_TRUE(c.bubble_map.find(s.node) != c.bubble_map.end()); // The bubble is registered
                res_entries.emplace_back(pos);
            }
            else if(is_site_exit) res_exits.emplace_back(pos);
        }
        catch(std::out_of_range) { // The node is not in any site; its site ID should be 0.
            EXPECT_EQ(site_ID, 0);
            continue;
        }
      }
    EXPECT_EQ(res_entries, expected_site_entry_points);
    EXPECT_EQ(res_exits, expected_site_exit_points);
}


class cov_G_Builder_nested_adjMarkers : public ::testing::Test {
protected:
    void SetUp() {
        // A nested string with adjacent variant markers
        // Namely due to: i)direct deletion and ii)double entry
        std::string prg_string{"[A,]A[[G,A]A,C]"};
        marker_vec v = prg_string_to_ints(prg_string);
        PRG_String p{v};
        p.process();
        c = cov_Graph_Builder{p};
        c.run();
    }
    cov_Graph_Builder c;
};

TEST_F(cov_G_Builder_nested_adjMarkers, adjMarkerWiring){
    //"[A,]A[[G,A]A,C]"
    covG_ptr entry;
   entry = c.bubble_starts.at(5);
   EXPECT_EQ(entry, c.random_access[0].node); // Consistent site numbering, sanity check
   auto& expected_exit = c.bubble_ends.at(5);
   // Expect direct edge between the site starting at index 0 and its site end
   EXPECT_EQ(entry->get_edges().back(), expected_exit);

    entry = c.bubble_starts.at(7);
    EXPECT_EQ(entry, c.random_access[5].node); // Consistent site numbering, sanity check
    auto& expected_next_entry = c.bubble_starts.at(9);
    // Expect direct edge between the site starting at index 5 and the site starting at index 6
    EXPECT_EQ(entry->get_edges()[0], expected_next_entry);
};
