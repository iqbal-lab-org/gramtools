/** @file
 * Test Level Genotyping (LG).
 * These are high-level tests:
 *  build a coverage graph & gram index, map reads to it, call a genotyper; all those are required to work.
 */

#include "gtest/gtest.h"

#include "genotype/quasimap/quasimap.hpp"
#include "genotype/infer/level_genotyping/runner.hpp"
#include "genotype/infer/output_specs/make_json.hpp"

#include "tests/genotype/infer/mocks.hpp"
#include "tests/test_resources/test_resources.hpp"

TEST(LevelGenotyping, Given2SiteNonNestedPRG_CorrectGenotypes){
    std::string prg{"AATAA5C6G6AA7C8G8AA"};
    Sequences kmers{encode_dna_bases("AA")};
    prg_setup setup;
    setup.setup_numbered_prg(prg, kmers);

    allele_vector gt_alleles;
    GenomicRead_vector reads;
    // Multiple reads going through 5:1 and 7:1
    for (int i = 0; i < 5; i++) reads.push_back(GenomicRead("Read", "AATAACAACAA", "???????????"));

    // One read going through 5:2 and 7:1
    reads.push_back(GenomicRead("ErrorRead", "AATAAGAACAA", "???????????"));

    setup.quasimap_reads(reads);

    LevelGenotyper genotyper(setup.prg_info.coverage_graph, setup.coverage.grouped_allele_counts,
                             setup.read_stats, Ploidy::Haploid);
    auto gt_recs = genotyper.get_genotyped_records();

    gt_alleles = gt_recs.at(siteID_to_index(5))->get_unique_genotyped_alleles();
    allele_vector expected_alleles{ Allele{"C", {5}, 0} };
    EXPECT_EQ(gt_alleles, expected_alleles);

    gt_alleles = gt_recs.at(0)->get_unique_genotyped_alleles();
    expected_alleles = allele_vector{ Allele{"C", {5}, 0} };
    EXPECT_EQ(gt_alleles, expected_alleles);
}

TEST(LevelGenotyping, Given2SiteNestedPRG_CorrectGenotypes){
    std::string prg{"AATAA[CCC[A,G],T]AA"};
    Sequences kmers{encode_dna_bases("AA")};
    prg_setup setup;
    setup.setup_bracketed_prg(prg, kmers);

    allele_vector gt_alleles;
    GenomicRead_vector reads;
    // Multiple reads going through first allele of each site
    for (int i = 0; i < 5; i++) reads.push_back(GenomicRead("Read", "AATAACCCGAA", "???????????"));
    // One read going through second allele of site 1 and first allele of site 2
    reads.push_back(GenomicRead("ErrorRead", "AATAATAA", "????????"));

    setup.quasimap_reads(reads);

    LevelGenotyper genotyper(setup.prg_info.coverage_graph, setup.coverage.grouped_allele_counts,
                             setup.read_stats, Ploidy::Haploid);
    auto gt_recs = genotyper.get_genotyped_records();

    gt_alleles = gt_recs.at(1)->get_unique_genotyped_alleles();
    allele_vector expected_alleles{
            Allele{"G", {5}, 1}
    };
    EXPECT_EQ(gt_alleles, expected_alleles);

    gt_alleles = gt_recs.at(0)->get_unique_genotyped_alleles();
    expected_alleles = allele_vector{
            Allele{"CCCG", {5, 5, 5, 5}, 0}
    };
    EXPECT_EQ(gt_alleles, expected_alleles);
}

TEST(LevelGenotyper, GivenPRGWithDirectDeletion_CorrectlyCalledEmptyAllele){
    std::string prg{"GGGGG[CCC,]GG"};
    Sequences kmers{encode_dna_bases("GG")};
    prg_setup setup;
    setup.setup_bracketed_prg(prg, kmers);

    allele_vector gt_alleles;
    GenomicRead_vector reads;
    // Reads going through direct deletion
    for (int i = 0; i < 5; i++) reads.push_back(GenomicRead("Read", "GGGGGG", "??????"));
    setup.quasimap_reads(reads);

    LevelGenotyper genotyper(setup.prg_info.coverage_graph, setup.coverage.grouped_allele_counts,
                             setup.read_stats, Ploidy::Haploid);
    auto gt_recs = genotyper.get_genotyped_records();

    gt_alleles = gt_recs.at(0)->get_unique_genotyped_alleles();
    allele_vector expected_alleles{
            Allele{"", {}, 1}
    };
    EXPECT_EQ(gt_alleles, expected_alleles);
}

class LG_SnpsNestedInTwoHaplotypes : public ::testing::Test {
protected:
    void SetUp(){
        std::string _prg{"ATCGGC[TC[A,G]TC,GG[T,G]GG]AT"};
        auto all_kmers = gram::generate_all_kmers(2);
        Sequences all_kmers_vector{all_kmers.begin(), all_kmers.end()};
        setup.setup_bracketed_prg(_prg, all_kmers_vector);

        // This read goes through 5:1 and 7:2
        for (int num_mapped{0}; num_mapped < 7; num_mapped++) reads.push_back(
                    GenomicRead("Read1", "ATCGGCTCGTCAT", ".............") );

        // This read goes through 5:2 and 9:2
        reads.push_back( GenomicRead("Read2", "ATCGGCGGG", ".........") );
    }

    void MapReadsAndHaploidGenotype(){
        setup.quasimap_reads(reads);
        LevelGenotyper genotyper(setup.prg_info.coverage_graph, setup.coverage.grouped_allele_counts,
                                 setup.read_stats, Ploidy::Haploid);
        gt_recs = genotyper.get_genotyped_records();
    }

    GenomicRead_vector reads;
    allele_vector gt_alleles;
    allele_vector expected_alleles;
    gt_sites gt_recs;
    prg_setup setup;
};

TEST_F(LG_SnpsNestedInTwoHaplotypes, MapNoReads_AllGenotypesAreNull){
    LevelGenotyper genotyper(setup.prg_info.coverage_graph, setup.coverage.grouped_allele_counts,
                             setup.read_stats, Ploidy::Haploid);
    gt_recs = genotyper.get_genotyped_records();

    for (auto const& gt_rec : gt_recs){
        EXPECT_TRUE(gt_rec->is_null());
    }
}

TEST_F(LG_SnpsNestedInTwoHaplotypes, MapReads_CorrectlyGenotypedSites){
    MapReadsAndHaploidGenotype();

    gt_alleles = gt_recs.at(siteID_to_index(5))->get_unique_genotyped_alleles();
    expected_alleles = allele_vector{
            Allele{"TCGTC", {7, 7, 7, 7, 7}, 0}
    };
    EXPECT_EQ(gt_alleles, expected_alleles);

    gt_alleles = gt_recs.at(siteID_to_index(7))->get_unique_genotyped_alleles();
    expected_alleles = allele_vector{
            Allele{"G", {7}, 1}
    };
    EXPECT_EQ(gt_alleles, expected_alleles);
}

TEST_F(LG_SnpsNestedInTwoHaplotypes, MapReads_CorrectlyInvalidatedSites){
    // Since we called 5:1, we should invalidate whatever lives on 5:2; which is site ID 9.
    MapReadsAndHaploidGenotype();

    EXPECT_TRUE(gt_recs.at(siteID_to_index(9))->is_null());

    auto site_result = gt_recs.at(siteID_to_index(9));
    auto json_result = make_json_site(site_result)->get_site();
    EXPECT_FLOAT_EQ(json_result.at("GT_CONF").at(0), 0.);
}

TEST(LevelGenotyperInvalidation, GivenChildMapAndCandidateHaplos_CorrectHaplosWithSites){
    // site 7 lives on haplogroup 0 of site 5, and sites 9 and 11 live on its haplogroup 1.
    parental_map par_map{
            {7, VariantLocus{5, FIRST_ALLELE}},
            {9, VariantLocus{5, FIRST_ALLELE + 1}},
            {11, VariantLocus{5, FIRST_ALLELE + 1}},
    };
    child_map child_m = build_child_map(par_map);
    LevelGenotyper g{child_m, gt_sites{}};

    AlleleIds expected_haplogroups{0, 1}; // Expected in 0-based
    auto haplos_with_sites = g.get_haplogroups_with_sites(5, {0, 1, 2, 3});
    EXPECT_EQ(haplos_with_sites, expected_haplogroups);

    auto empty_query = g.get_haplogroups_with_sites(7, {0, 1, 2, 3});
    EXPECT_EQ(empty_query, AlleleIds{});
}

using ::testing::InSequence;
using ::testing::Return;
TEST(LevelGenotyperInvalidation, GivenNestingStructure_CorrectGenotypeNullifying){
    parental_map par_map{
            {7, VariantLocus{5, FIRST_ALLELE}},
            {9, VariantLocus{7, FIRST_ALLELE + 1}},
    };
    child_map child_m = build_child_map(par_map);

    gt_sites sites(3);
    auto site1 = std::make_shared<MockGenotypedSite>();
    site1->set_num_haplogroups(5);
    sites.at(1) = site1;

    // SiteID 9 will get nulled by site 7.
    // Then when site 5 nulls site 7, I want site 9 to signal it is already nulled.
    auto site2 = std::make_shared<MockGenotypedSite>();
    site2->set_num_haplogroups(5);
    sites.at(2) = site2;

    LevelGenotyper g{child_m, sites};

    EXPECT_FALSE(site2->is_null());
    // I want site 9 to be invalidated in this call
    g.invalidate_if_needed(7, AlleleIds{1});
    EXPECT_TRUE(site2->is_null());

    EXPECT_FALSE(site1->is_null());
    // And this call to invalidate site 7, without attempting to invalidate site 9 which was already so by above call.
    g.invalidate_if_needed(5, AlleleIds{0});
    EXPECT_TRUE(site1->is_null());
}
