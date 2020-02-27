#include "genotype/infer/personalised_reference.hpp"
#include "prg/coverage_graph.hpp"
#include "gtest/gtest.h"
#include "mocks.hpp"
#include "tests/common/common.hpp"

using namespace gram::genotype;

class Alleles_To_Paste : public ::testing::Test{
protected:
    void SetUp(){
       site->set_alleles(all_alleles);
    }

    gt_site_ptr site = std::make_shared<MockGenotypedSite>();

    allele_vector all_alleles{
            Allele{"ATA", {0, 0, 0}, 0},
            Allele{"TTA", {0, 0, 0}, 1},
            Allele{"TTT", {0, 0, 0}, 2},
    };
};

TEST_F(Alleles_To_Paste, GivenInconsistentPloidy_Throws){
   site->set_genotype(GtypedIndices{0, 1});
   EXPECT_THROW(get_all_alleles_to_paste(site, 3), InconsistentPloidyException);
}

TEST_F(Alleles_To_Paste, GivenGtype_CorrectAlleles){
    site->set_genotype(GtypedIndices{0, 2});
    auto res = get_all_alleles_to_paste(site, 2);
    allele_vector expected{all_alleles.at(0), all_alleles.at(2)};
    EXPECT_EQ(res, expected);
}

TEST_F(Alleles_To_Paste, GivenNullGtype_CorrectAlleles){
    site->set_genotype(false);
    auto res = get_all_alleles_to_paste(site, 3);
    allele_vector expected{all_alleles.at(0), all_alleles.at(0), all_alleles.at(0)};
    EXPECT_EQ(res, expected);
}

class Personalised_Reference : public ::testing::Test{
protected:
    void SetUp(){
        std::string linear_prg{"AT[CG[C,G]T,C]TT[AT,TT]"};
        PRG_String prg{prg_string_to_ints(linear_prg)};
        g = coverage_Graph{prg};
        graph_root = g.root;
        covG_ptrPair bubble;

        auto site1 = std::make_shared<MockGenotypedSite>();
        site1->set_alleles(allele_vector{
            Allele{"CGCT", {}, 0},
            Allele{"CGGT", {}, 0},
            Allele{"C", {}, 1},
        });
        bubble = get_bubble_nodes(g.bubble_map, 5);
        site1->set_site_end_node(bubble.second);

        // This, being nested, should get systematically skipped
        auto site2 = std::make_shared<MockGenotypedSite>();
        site2->set_alleles(allele_vector{
                Allele{"C", {}}, Allele{"G", {}},
        });
        bubble = get_bubble_nodes(g.bubble_map, 7);
        site2->set_site_end_node(bubble.second);

        auto site3 = std::make_shared<MockGenotypedSite>();
        site3->set_alleles(allele_vector{
                Allele{"AT", {}},
                Allele{"TT", {}},
        });
        bubble = get_bubble_nodes(g.bubble_map, 9);
        site3->set_site_end_node(bubble.second);

        sites = gt_sites{site1, site2, site3};
    }
    coverage_Graph g;
    covG_ptr graph_root;
    gt_sites sites;
};

TEST_F(Personalised_Reference, GivenAllNullGts_CorrectInferredRef){
    // When all gts are null, ploidy is set to 1
    for (auto const& site : sites) site->set_genotype(false);
    auto result_map = get_personalised_ref(graph_root, sites);
    auto result = *result_map.begin();
    std::string expected{"ATCGCTTTAT"};
    EXPECT_EQ(result.get_sequence(), expected);
}

TEST_F(Personalised_Reference, GivenHaploidGts_CorrectInferredRef){
    sites.at(0)->set_genotype(GtypedIndices{2});
    sites.at(2)->set_genotype(GtypedIndices{1});
    auto result_map = get_personalised_ref(graph_root, sites);
    auto result = *result_map.begin();
    std::string expected{"ATCTTTT"};
    EXPECT_EQ(result.get_sequence(), expected);
}

TEST_F(Personalised_Reference, GivenHetDiploidGts_CorrectTwoInferredRefs){
    sites.at(0)->set_genotype(GtypedIndices{1, 2});
    sites.at(2)->set_genotype(GtypedIndices{0, 1});
    auto result_map = get_personalised_ref(graph_root, sites);
    EXPECT_EQ(result_map.size(), 2);
    auto iterator = result_map.begin();

    auto first_ref = *iterator;
    std::string expected_1{"ATCGGTTTAT"};
    EXPECT_EQ(first_ref.get_sequence(), expected_1);

    iterator++;
    auto second_ref = *iterator;
    std::string expected_2{"ATCTTTT"};
    EXPECT_EQ(second_ref.get_sequence(), expected_2);
}


TEST_F(Personalised_Reference, GivenHetSameGts_CorrectSingleInferredRef) {
    sites.at(0)->set_genotype(GtypedIndices{0, 0});
    sites.at(2)->set_genotype(GtypedIndices{1, 1});
    auto result_vec = get_personalised_ref(graph_root, sites);
    EXPECT_EQ(result_vec.size(), 2);

    unique_Fastas result_map{result_vec.begin(), result_vec.end()};

    auto first_ref = *result_map.begin();
    std::string expected_1{"ATCGCTTTTT"};
    EXPECT_EQ(first_ref.get_sequence(), expected_1);
}
