#include <filesystem>
#include "gtest/gtest.h"

#include "prg/coverage_graph.hpp"
#include "../test_resources/test_resources.hpp"

using namespace gram;
namespace fs = std::filesystem;

auto const test_data_dir = fs::path(__FILE__).parent_path().parent_path() / "test_data";

auto first{FIRST_ALLELE};
auto unkn{ALLELE_UNKNOWN};

/*
 * -----------------------
 * `PRG String` tests
 * -----------------------
 */

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

/**
 * Given inconsistent PRGs, prg_string/cov_graph construction fails
 */

TEST(InconsistentPRG_duplicate_site_markers, fails){
    marker_vec t{5,1,6,2,6,2,5,1,6,3,6}; // "[A,C]C[A,G]"
    EXPECT_THROW(PRG_String{t}, std::runtime_error);
}

TEST(InconsistentPRG_site_with_no_alleles, fails){
    marker_vec m{5,6,2,7,1,8,3,8}; // "[]C[A,G]"
    PRG_String s{m};
    EXPECT_THROW(cov_Graph_Builder{s}, std::runtime_error);
}

TEST(InconsistentPRG_site_one_allele, fails) {
    marker_vec m{5,2,6,2,7,1,8,3,8}; // "[C]C[A,G]"
    PRG_String s{m};
    EXPECT_THROW(cov_Graph_Builder{s}, std::runtime_error);
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
    // Expecting to find a single entry, for the single nested site, pointing to siteID 5 & alleleID 3.
    parental_map expected {
            {7 , VariantLocus{5, 2} }
    };
    EXPECT_EQ(c.par_map, expected);
}

// Test that the node site & allele IDs are correct
TEST_F(cov_G_Builder_nested, SiteAndAllele_IDs){
    //"[A,AA,A[A,C]A]C[AC,C]G"
    auto const& rand_access = c.random_access;
    std::vector<VariantLocus> expected{
            {5, unkn}, {5, first}, {5, unkn}, {5, first + 1}, {5, first + 1},
            {5, unkn}, {5, first + 2}, {7, unkn}, {7, first}, {7, unkn},
            {7, first + 1}, {7, unkn}, {5, first + 2}, {5, unkn},
            {0, unkn}, {9, unkn}, {9, first}, {9, first},
            {9, unkn}, {9,first + 1}, {9, unkn}, {0, unkn}
    };
    std::vector<VariantLocus> res{rand_access.size(), VariantLocus()};
    int pos = 0;
    for (auto const&s : rand_access){
       res[pos].first = s.node->get_site_ID();
        res[pos].second = s.node->get_allele_ID();
        pos++;
    }

    EXPECT_EQ(res, expected);
}

// Test that the size of the nodes is correct
TEST_F(cov_G_Builder_nested, NodeSizes) {
    //"[A,AA,A[A,C]A]C[AC,C]G"
    auto const &rand_access = c.random_access;
    // This test queries UNIQUE nodes, so we will skip "," which point to bubble start node,
    // and sequence continuation for nodes with size > 1
    std::vector<int> expected{
        0, 1, 2, 1, 0, 1, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1
    };
    std::vector<int> res(expected.size(), -1);
    std::set<Marker> seen_entries; // For skipping bubble entry nodes
    int pos = 0;
    covG_ptr prev = nullptr; // For skipping consecutive nucleotides
    for (auto const&s : rand_access) {
        // Test for skipping site entry points
        if (c.bubble_map.find(s.node) != c.bubble_map.end()){
            if (seen_entries.find(s.node->get_site_ID()) != seen_entries.end()) continue;
            else seen_entries.insert(s.node->get_site_ID());
        }
        if (s.node == prev) continue;
        auto sequence_size = s.node->get_sequence_size();

        // Test there is as much allocated per base coverage as there are characters in the sequence node
        // if we are in variant site. Outside variant sites we do not genotype so do not allocate/record coverage.
        if (s.node->is_in_bubble()) EXPECT_EQ(s.node->get_coverage_space(), sequence_size);

        res[pos++] = sequence_size;
        prev = s.node;
    }
    EXPECT_EQ(res,expected);
}

// Test that the node positions are correct
TEST_F(cov_G_Builder_nested, SequencePositions) {
    //"[A,AA,A[A,C]A]C[AC,C]G"
    auto const &rand_access = c.random_access;
    // There is one position per index in the PRG string
    // The positions continually refer to the FIRST allele in each site, 0-based.
    std::vector<std::size_t> expected{
        0, 0, 0, 0, 0, 0, 0, // First site here
            1, 1, 1, 1, 2, // Second site FULL here
        2, 1,  // First site END here
        1, // Invariant 'C'
        2, 2, 2, 2, 2, 4, // Third site FULL here
        4 // Invariant 'G'
        };
    std::vector<std::size_t> res(expected.size(), 0);
    int pos = 0;
    for (auto const &s : rand_access) {
        res[pos++] = s.node->get_pos();
    }
    EXPECT_EQ(res, expected);
}


// Test that bubble entry and exit points are correctly identified
TEST_F(cov_G_Builder_nested, Bubble_Positions){
    //"[A,AA,A[A,C]A]C[AC,C]G
    using i_v = std::vector<int>; // Stores indexes into PRG string
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
        Marker site_ID = s.node->get_site_ID();
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
        std::string prg_string{"[A,]A[[G,A]A,C,T]"};
        marker_vec v = prg_string_to_ints(prg_string);
        PRG_String p{v};
        c = cov_Graph_Builder{p};
    }
    cov_Graph_Builder c;
};

TEST_F(cov_G_Builder_nested_adjMarkers, adjMarkerWiring){
    //"[A,]A[[G,A]A,C,T]"
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

TEST_F(cov_G_Builder_nested_adjMarkers, bubbleOrdering){
    //"[A,]A[[G,A]A,C,T]"
    /*
     * Tests the bubbles are in the right order.
     * Because of the double entry, two bubbles have the same POS,
     * so make sure the more nested (child) bubble occurs before its parent.
     * This is needed for nested genotyping.
     */
    std::vector<std::size_t> expected{2, 1, 0};
    std::vector<std::size_t> site_indices;
    for (auto const& bubble_pair: c.bubble_map)
        site_indices.push_back(siteID_to_index(bubble_pair.first->get_site_ID()));
    EXPECT_EQ(site_indices, expected);
}

// Tests the target mapping is correct
TEST_F(cov_G_Builder_nested_adjMarkers, targetEntries){
    //"[A,]A[[G,A]A,C,T]"
    /**
     * First, check that nucleotide positions just after a marker target the site and allele markers
     */
    marker_vec expected_site_targets = {
            0, 5, 0, 0, 6, 0, 0, 9, 0, 10, 0, 10, 0, 8, 0, 8, 0
    };
    AlleleIds expected_allele_targets = {
            unkn, first, unkn, unkn, unkn, unkn, unkn,
            first, unkn, first + 1, unkn, first, unkn, first + 1, unkn,
            first + 2, unkn
    };

    marker_vec site_results(expected_site_targets.size(), 0);
    AlleleIds allele_results(expected_site_targets.size(), unkn);
    int pos = -1;
    for (auto& e : c.random_access){
        pos++;
        site_results[pos] = e.target.first;
        allele_results[pos] = e.target.second;
    }
    EXPECT_EQ(site_results, expected_site_targets);
    EXPECT_EQ(allele_results, expected_allele_targets);

    /**
     * Second, check that adjacent variant markers get correct entries in the target map
     */
    std::vector<targeted_marker> v;
    Marker seed;
    target_m expected_map;
    // First add in the direct deletion at pos 3
    seed = 6;
    v.emplace_back(targeted_marker{5, first + 1});
    expected_map.insert(std::make_pair(seed, v));
    v.clear();

    // Second add the double start at pos 6
    seed = 9;
    v.emplace_back(targeted_marker{7,ALLELE_UNKNOWN});
    expected_map.insert(std::make_pair(seed, v));
    v.clear();

    EXPECT_EQ(c.target_map, expected_map);
}

// Test for the number of sites, and that each "," character amounts to returning to the site entry point
TEST_F(cov_G_Builder_nested_adjMarkers, num_Bubbles) {
    //"[A,]A[[G,A]A,C,T]"

    // Will record how many times each site entry node has been traversed
    std::unordered_map<Marker, int> seen_entries;
    std::unordered_map<Marker, int> expected{
            {5, 1},
            {7, 2},
            {9, 1}
    };

    for (auto const &s : c.random_access) {
        if (c.bubble_map.find(s.node) != c.bubble_map.end()) {
            if (seen_entries.find(s.node->get_site_ID()) != seen_entries.end()) {
                // Case: at the entry node for the at least second time
                seen_entries.at(s.node->get_site_ID())++;
            }
            // Case: at the entry node for the first time
            else seen_entries.insert({s.node->get_site_ID(), 0});
        }
    }

    EXPECT_EQ(seen_entries, expected);
}

// Test that the parental map deals with adjacent markers
TEST_F(cov_G_Builder_nested_adjMarkers, ParentalMap) {
    //"[A,]A[[G,A]A,C,T]"
    parental_map expected {
            {9 , VariantLocus{7, first} }
    };
    EXPECT_EQ(c.par_map, expected);
}

TEST(coverage_Graph, Nestedness){
    std::string prg{"ATCG[GC,G]A[AT,T]A"};
    marker_vec v = prg_string_to_ints(prg);
    PRG_String p{v};
    coverage_Graph g{p};

    EXPECT_FALSE(g.is_nested);

    std::string nested_prg{"[A,]A[[G,A]A,C,T]"};
    marker_vec nested_v = prg_string_to_ints(nested_prg);
    PRG_String nested_p{nested_v};
    coverage_Graph nested_g{nested_p};

    EXPECT_TRUE(nested_g.is_nested);
}

TEST(coverage_Graph, SequencePositions) {
    // Check POS is based on first (REF) allele of each site
    std::string prg{"ATCG[G[A,CCC]C,G]A[AT,T]A"};
    marker_vec v = prg_string_to_ints(prg);
    PRG_String p{v};
    coverage_Graph g{p};

    auto bubble_5 = get_bubble_nodes(g.bubble_map, 5);
    EXPECT_EQ(4, bubble_5.first->get_pos());

    auto bubble_7 = get_bubble_nodes(g.bubble_map, 7);
    EXPECT_EQ(5, bubble_7.first->get_pos());

    auto bubble_9 = get_bubble_nodes(g.bubble_map, 9);
    EXPECT_EQ(8, bubble_9.first->get_pos());
}

TEST(coverage_Graph, SequencePositions2) {
    // Check POS is updated for first allele only in sites with nesting
    std::string prg{"ATCG[G[A,CCC]C,GGG[AAA,C]]AA[T,C]"};
    marker_vec v = prg_string_to_ints(prg);
    PRG_String p{v};
    coverage_Graph g{p};

    auto bubble_5 = get_bubble_nodes(g.bubble_map, 5);
    EXPECT_EQ(4, bubble_5.first->get_pos());

    auto bubble_7 = get_bubble_nodes(g.bubble_map, 7);
    EXPECT_EQ(5, bubble_7.first->get_pos());

    auto bubble_9 = get_bubble_nodes(g.bubble_map, 9);
    EXPECT_EQ(7, bubble_9.first->get_pos());

    auto bubble_11 = get_bubble_nodes(g.bubble_map, 11);
    EXPECT_EQ(9, bubble_11.first->get_pos());
}


// Make a coverage graph, serialise it to disk, reload into another coverage graph,
// and test the two are equal (provided equality has been properly defined).
TEST(coverage_Graph, Serialisation){
    std::string prg_string{"[A,]A[[G,A]A,C,T]"};
    marker_vec v = prg_string_to_ints(prg_string);
    PRG_String p{v};
    coverage_Graph serialised_cov_G{p};

    fs::path path(test_data_dir / "tmp.ar");

    // Dump to disk
    {
        std::ofstream ofs{path.generic_string()};
        EXPECT_TRUE(ofs); // Can make this file
        boost::archive::binary_oarchive oa{ofs};
        oa << serialised_cov_G;
    }
    EXPECT_TRUE(fs::exists(path)); // Have made this file

    // Load from disk
    coverage_Graph loaded_cov_G;
    std::ifstream ifs{path.generic_string()};
    boost::archive::binary_iarchive ia{ifs};
    ia >> loaded_cov_G;

    EXPECT_TRUE(serialised_cov_G == loaded_cov_G);
}

TEST(Target_map, EvenIsEntry_OddIsExit){
    std::string prg_string{"[A,[A,C,G]A,C]"};
    marker_vec v = prg_string_to_ints(prg_string);
    PRG_String p{v};
    auto c = cov_Graph_Builder{p};

    std::vector<targeted_marker> targets{
        targeted_marker{5, ALLELE_UNKNOWN}
    };
    Marker seed{7};

    target_m expected_map;
    expected_map.insert(std::make_pair(seed, targets));

    EXPECT_EQ(c.target_map, expected_map);
}
