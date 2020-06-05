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


class TestInduceGenotypes_ThreadSequence : public ::testing::Test{
protected:
    coverage_Graph g;
    covG_ptr graph_end;
    TestInduceGenotypes_ThreadSequence(){
        auto encoded_prg = prg_string_to_ints("AA[A,C,G]TG[AC,[G,T]CA]CCC");
        PRG_String p{encoded_prg};
        g = coverage_Graph{p};
        auto cur_Node = g.root;
        while (cur_Node->get_num_edges() > 0)
            cur_Node = cur_Node->get_edges().at(0);
        graph_end = cur_Node;
    }
};

TEST_F(TestInduceGenotypes_ThreadSequence, GivenSequenceNotInGraph_ThrowsError){
    std::string absent_sequence{"AACTGACTTT"};
    EXPECT_THROW(thread_sequence(g.root, absent_sequence), NoEndpoints);
}

TEST_F(TestInduceGenotypes_ThreadSequence, GivenSequenceInGraphButIncomplete_ThrowsError){
    std::string sequence{"AACTGACC"};
    EXPECT_THROW(thread_sequence(g.root, sequence), NoEndpoints);
}

TEST_F(TestInduceGenotypes_ThreadSequence, GivenSequenceInGraphAndComplete_GetSingleEndpoint){
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

TEST(InduceGenotypes_MakeNullSites, SitesAreNullGTAndHaveRefSeqOnly){
    auto encoded_prg = prg_string_to_ints("AT[C,C[A,T]]GG");
    auto g = coverage_Graph{PRG_String{encoded_prg}};
    auto sites = make_nulled_sites(g);
    for (auto const& site : sites){
        EXPECT_TRUE(site->is_null());
        EXPECT_EQ(site->get_alleles().size(), 1);
    }
    auto site1 = sites.at(0);
    EXPECT_EQ(site1->get_alleles().at(0).sequence, "C");

    auto site2 = sites.at(1);
    EXPECT_EQ(site2->get_alleles().at(0).sequence, "A");
}


class TestInduceGenotypes_ApplyGenotypes : public ::testing::Test{
protected:
    gt_sites sites;
    coverage_Graph g;
    void SetUp(){
        auto encoded_prg = prg_string_to_ints("AT[,C,GG]AA[TA,AA,G[GG,GGG]A,]CA");
        g = coverage_Graph{PRG_String{encoded_prg}};
        sites = make_nulled_sites(g);
    }
};

TEST_F(TestInduceGenotypes_ApplyGenotypes, GivenRefThreadedSeq_CorrectGenotypedSites){
    auto endpoint = thread_sequence(g.root, "ATAATACA");
    apply_genotypes(endpoint, sites);

    for (auto const& site : gt_sites{sites.begin(), sites.begin() + 2}){
        EXPECT_FALSE(site->is_null());
        auto result = site->get_all_gtype_info();
        EXPECT_EQ(result.alleles.size(), 1);
        EXPECT_EQ(result.genotype, GtypedIndices{0});
        EXPECT_EQ(result.haplogroups, AlleleIds{0});
    }

    auto alleles = sites.at(0)->get_alleles();
    EXPECT_EQ(alleles.at(0).sequence, "");

    alleles = sites.at(1)->get_alleles();
    EXPECT_EQ(alleles.at(0).sequence, "TA");

    EXPECT_TRUE(sites.at(2)->is_null());
}

TEST_F(TestInduceGenotypes_ApplyGenotypes, GivenNonRefThreadedSeq_CorrectGenotypedSites) {
    auto endpoint = thread_sequence(g.root, "ATCAAGGGGACA");
    apply_genotypes(endpoint, sites);
    std::vector<std::string> expected_seqs;
    AlleleIds expected_ids;

    for (auto const& site : sites){
        EXPECT_FALSE(site->is_null());
        auto result = site->get_all_gtype_info();
        EXPECT_EQ(result.alleles.size(), 2);
        EXPECT_EQ(result.genotype, GtypedIndices{1});
        EXPECT_EQ(result.haplogroups.size(), 2);

        expected_seqs.push_back(result.alleles.back().sequence);
        expected_ids.push_back(result.haplogroups.back());
    }

    EXPECT_EQ(expected_seqs, std::vector<std::string>({"C", "GGGGA", "GGG"}));
    EXPECT_EQ(expected_ids, AlleleIds({1, 2, 1}));
}
