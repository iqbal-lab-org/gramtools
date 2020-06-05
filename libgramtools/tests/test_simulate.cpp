#include "gtest/gtest.h"

#include "simulate/simulate.hpp"
#include "simulate/induce_genotypes.hpp"
#include "prg/linearised_prg.hpp"
#include "prg/coverage_graph.hpp"

#include "test_resources/mocks.hpp"

using namespace gram::simulate;

class MakeRandomGenotypedSite : public ::testing::Test{
protected:
    allele_vector alleles{
            Allele{"CTCGG", {}},
            Allele{"CG", {}},
            Allele{"CT", {}},
    };
    MockRandomGenerator rand;
};

using ::testing::Return;

TEST_F(MakeRandomGenotypedSite, GivenPickZerothAllele_CorrectSite){
   EXPECT_CALL(rand, generate(0, 2))
   .WillOnce(Return(0));
   auto site = make_randomly_genotyped_site(&rand, alleles, true);
   allele_vector expected_als{alleles.begin(), alleles.begin() + 1};
   EXPECT_EQ(site->get_alleles(), expected_als);

   GtypedIndices expected_gts{0};
   EXPECT_EQ(site->get_genotype(), expected_gts);

   EXPECT_EQ(site->get_num_haplogroups(), 3);
}

TEST_F(MakeRandomGenotypedSite, GivenPickSecondAllele_CorrectSite){
    EXPECT_CALL(rand, generate(0, 2))
            .WillOnce(Return(2));
    auto site = make_randomly_genotyped_site(&rand, alleles, true);
    allele_vector expected_als{
        alleles.at(0),
        alleles.at(2)
    };
    EXPECT_EQ(site->get_alleles(), expected_als);

    GtypedIndices expected_gts{1}; // Is rescaled
    EXPECT_EQ(site->get_genotype(), expected_gts);
}

TEST_F(MakeRandomGenotypedSite, GivenIgnoreREFAllele_CorrectSite){
    EXPECT_CALL(rand, generate(1, 2))
            .WillOnce(Return(1));
    auto site = make_randomly_genotyped_site(&rand, alleles, false); // false: do not consider REF
    allele_vector expected_als{
            alleles.at(0),
            alleles.at(1)
    };
    EXPECT_EQ(site->get_alleles(), expected_als);
}

class TestInduceGenotypes : public ::testing::Test{
protected:
    coverage_Graph g;
    covG_ptr graph_end;
    TestInduceGenotypes(){
        auto encoded_prg = prg_string_to_ints("AA[A,C,G]TG[AC,[G,T]CA]CCC");
        PRG_String p{encoded_prg};
        g = coverage_Graph{p};
        auto cur_Node = g.root;
        while (cur_Node->get_num_edges() > 0)
            cur_Node = cur_Node->get_edges().at(0);
        graph_end = cur_Node;
    }
};

TEST_F(TestInduceGenotypes, GivenSequenceNotInGraph_ThrowsError){
    std::string absent_sequence{"AACTGACTTT"};
    EXPECT_THROW(thread_sequence(g.root, absent_sequence), NoEndpoints);
}

TEST_F(TestInduceGenotypes, GivenSequenceInGraphButIncomplete_ThrowsError){
    std::string sequence{"AACTGACC"};
    EXPECT_THROW(thread_sequence(g.root, sequence), NoEndpoints);
}

TEST_F(TestInduceGenotypes, GivenSequenceInGraphAndComplete_GetSingleEndpoint){
    std::string goodseq1{"AACTGACCCC"};
    auto result = thread_sequence(g.root, goodseq1);
    EXPECT_EQ(result->get_prg_node(), graph_end);
    EXPECT_EQ(result->get_offset(), 10);

    std::string goodseq2{"AAATGGCACCC"};
    result = thread_sequence(g.root, goodseq2);
    EXPECT_EQ(result->get_prg_node(), graph_end);
    EXPECT_EQ(result->get_offset(), 11);
}

TEST(InduceGenotypesAmbiguous, GivenMultiplePossiblePathsInPrg_ThrowsError){
    /*
     * This highlights that the prg construction process needs to be sequence coherent:
     * Same sequences should not be produceable inside a site, or across full PRG.
     */
    auto ambiguous_prg = prg_string_to_ints("AA[A,AA]A[AA,A]T");
    // An equivalent, unambiguous prg is "AA[AA,AAA,AAAA]AT"
    auto g = coverage_Graph{PRG_String{ambiguous_prg}};
    EXPECT_THROW(thread_sequence(g.root, "AAAAAAT"),TooManyEndpoints);

    ambiguous_prg = prg_string_to_ints("AT[CA,C[C,A]]GG");
    g = coverage_Graph{PRG_String{ambiguous_prg}};
    EXPECT_THROW(thread_sequence(g.root, "ATCAGG"),TooManyEndpoints);
}
