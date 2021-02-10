/** @file
 * Test Level Genotyping (LG).
 * These are high-level tests:
 *  build a coverage graph & gram index, map reads to it, call a genotyper; all
 * those are required to work.
 */

#include "../../../test_resources/test_resources.hpp"
#include "../mocks.hpp"
#include "genotype/infer/level_genotyping/runner.hpp"
#include "genotype/infer/output_specs/make_json.hpp"
#include "gtest/gtest.h"

TEST(LevelGenotyping, Given2SiteNonNestedPRG_CorrectGenotypes) {
  std::string prg{"AATAA5C6G6AA7C8G8AA"};
  prg_setup setup;
  setup.setup_numbered_prg(prg);

  GenomicRead_vector reads;
  // Multiple reads going through 5:1 and 7:1
  for (int i = 0; i < 5; i++)
    reads.push_back(GenomicRead("Read", "AATAACAACAA", "???????????"));
  // One read going through 5:2 and 7:1
  reads.push_back(GenomicRead("ErrorRead", "AATAAGAACAA", "???????????"));

  setup.quasimap_reads(reads);

  LevelGenotyper genotyper(setup.prg_info.coverage_graph,
                           setup.coverage.grouped_allele_counts,
                           setup.read_stats, Ploidy::Haploid);
  auto gt_recs = genotyper.get_genotyped_records();

  allele_vector gt_alleles =
      gt_recs.at(siteID_to_index(5))->get_unique_genotyped_alleles();
  EXPECT_EQ(gt_alleles, (allele_vector{Allele{"C", {5}, 0}}));

  gt_alleles = gt_recs.at(0)->get_unique_genotyped_alleles();
  EXPECT_EQ(gt_alleles, (allele_vector{Allele{"C", {5}, 0}}));
}

TEST(LevelGenotyping, Given2SiteNestedPRG_CorrectGenotypes) {
  std::string prg{"AATAA[CCC[A,G],T]AA"};
  prg_setup setup;
  setup.setup_bracketed_prg(prg);

  GenomicRead_vector reads;
  // Multiple reads going through first allele of each site
  for (int i = 0; i < 5; i++)
    reads.push_back(GenomicRead("Read", "AATAACCCGAA", "???????????"));
  // One read going through second allele of site 1 and first allele of site 2
  reads.push_back(GenomicRead("ErrorRead", "AATAATAA", "????????"));

  setup.quasimap_reads(reads);

  LevelGenotyper genotyper(setup.prg_info.coverage_graph,
                           setup.coverage.grouped_allele_counts,
                           setup.read_stats, Ploidy::Haploid);
  auto gt_recs = genotyper.get_genotyped_records();

  allele_vector gt_alleles = gt_recs.at(1)->get_unique_genotyped_alleles();
  EXPECT_EQ(gt_alleles, (allele_vector{Allele{"G", {5}, 1}}));

  gt_alleles = gt_recs.at(0)->get_unique_genotyped_alleles();
  EXPECT_EQ(gt_alleles, (allele_vector{Allele{"CCCG", {5, 5, 5, 5}, 0}}));
}

TEST(LevelGenotyper, GivenPRGWithDirectDeletion_CorrectlyCalledEmptyAllele) {
  std::string prg{"GGGGG[CCC,]GG"};
  prg_setup setup;
  setup.setup_bracketed_prg(prg);

  allele_vector gt_alleles;
  GenomicRead_vector reads;
  // Reads going through direct deletion
  for (int i = 0; i < 5; i++)
    reads.push_back(GenomicRead("Read", "GGGGGG", "??????"));
  setup.quasimap_reads(reads);

  LevelGenotyper genotyper(setup.prg_info.coverage_graph,
                           setup.coverage.grouped_allele_counts,
                           setup.read_stats, Ploidy::Haploid);
  auto gt_recs = genotyper.get_genotyped_records();

  gt_alleles = gt_recs.at(0)->get_unique_genotyped_alleles();
  allele_vector expected_alleles{Allele{"", {}, 1}};
  EXPECT_EQ(gt_alleles, expected_alleles);
}

class LG_SnpsNestedInTwoHaplotypes : public ::testing::Test {
 protected:
  void SetUp() {
    std::string _prg{"ATCGGC[TC[A,G]TC,GG[T,G]GG]AT"};
    setup.setup_bracketed_prg(_prg);

    // This read goes through 5:0 and 7:1
    for (int num_mapped{0}; num_mapped < 7; num_mapped++)
      reads.push_back(GenomicRead("Read1", "ATCGGCTCGTCAT", "............."));

    // This read goes through 5:1 and 9:1
    reads.push_back(GenomicRead("Read2", "ATCGGCGGG", "........."));
  }

  void MapReadsAndHaploidGenotype() {
    setup.quasimap_reads(reads);
    LevelGenotyper genotyper(setup.prg_info.coverage_graph,
                             setup.coverage.grouped_allele_counts,
                             setup.read_stats, Ploidy::Haploid);
    gt_recs = genotyper.get_genotyped_records();
  }

  GenomicRead_vector reads;
  allele_vector gt_alleles;
  allele_vector expected_alleles;
  gt_sites gt_recs;
  prg_setup setup;
};

TEST_F(LG_SnpsNestedInTwoHaplotypes, MapNoReads_AllGenotypesAreNull) {
  LevelGenotyper genotyper(setup.prg_info.coverage_graph,
                           setup.coverage.grouped_allele_counts,
                           setup.read_stats, Ploidy::Haploid);
  gt_recs = genotyper.get_genotyped_records();

  for (auto const& gt_rec : gt_recs) {
    EXPECT_TRUE(gt_rec->is_null());
  }
}

TEST_F(LG_SnpsNestedInTwoHaplotypes, MapReads_CorrectlyGenotypedSites) {
  MapReadsAndHaploidGenotype();

  gt_alleles = gt_recs.at(siteID_to_index(5))->get_unique_genotyped_alleles();
  expected_alleles = allele_vector{Allele{"TCGTC", {7, 7, 7, 7, 7}, 0}};
  EXPECT_EQ(gt_alleles, expected_alleles);

  gt_alleles = gt_recs.at(siteID_to_index(7))->get_unique_genotyped_alleles();
  expected_alleles = allele_vector{Allele{"G", {7}, 1}};
  EXPECT_EQ(gt_alleles, expected_alleles);
}

TEST_F(LG_SnpsNestedInTwoHaplotypes, MapReads_CorrectlyInvalidatedSites) {
  // Since we called 5:1, we should invalidate whatever lives on 5:2; which is
  // site ID 9.
  MapReadsAndHaploidGenotype();

  EXPECT_TRUE(gt_recs.at(siteID_to_index(9))->is_null());

  auto site_result = gt_recs.at(siteID_to_index(9));
  auto json_result = make_json_site(site_result)->get_site();
  EXPECT_FLOAT_EQ(json_result.at("GT_CONF").at(0), 0.);
}

TEST(GCPSimulation, GivenDifferentNumGenotypedSites_ConsistentNumConfidences) {
  auto l_stats = LevelGenotyper::make_l_stats(20, 10, 0.1);
  Ploidy ploidy{Ploidy::Haploid};

  gt_sites sites(CONF_DISTRIB_SIZE);
  for (int i{0}; i < CONF_DISTRIB_SIZE; i++) {
    auto site = std::make_shared<LevelGenotypedSite>();
    site->set_gt_conf(10);
    sites.at(i) = std::static_pointer_cast<gt_site>(site);
  }
  auto confidences = LevelGenotyper::get_gtconf_distrib(sites, l_stats, ploidy);
  EXPECT_EQ(CONF_DISTRIB_SIZE, confidences.size());
  std::set<double> unique(confidences.begin(), confidences.end());
  EXPECT_EQ(1, unique.size());

  sites.resize(10);
  confidences = LevelGenotyper::get_gtconf_distrib(sites, l_stats, ploidy);
  EXPECT_EQ(CONF_DISTRIB_SIZE, confidences.size());
}

TEST(LevelGenotyperInvalidation,
     GivenChildMapAndCandidateHaplos_CorrectHaplosWithSites) {
  // site 7 lives on haplogroup 0 of site 5, and sites 9 and 11 live on its
  // haplogroup 1.
  parental_map par_map{
      {7, VariantLocus{5, FIRST_ALLELE}},
      {9, VariantLocus{5, FIRST_ALLELE + 1}},
      {11, VariantLocus{5, FIRST_ALLELE + 1}},
  };
  child_map child_m = build_child_map(par_map);
  LevelGenotyper g{child_m, gt_sites{}};

  AlleleIds expected_haplogroups{0, 1};  // Expected in 0-based
  auto haplos_with_sites = g.get_haplogroups_with_sites(5, {0, 1, 2, 3});
  EXPECT_EQ(haplos_with_sites, expected_haplogroups);

  auto empty_query = g.get_haplogroups_with_sites(7, {0, 1, 2, 3});
  EXPECT_EQ(empty_query, AlleleIds{});
}

class LevelGenotyperPropagation : public ::testing::Test {
 protected:
  parental_map par_map{
      {7, VariantLocus{5, FIRST_ALLELE}},
      {9, VariantLocus{7, FIRST_ALLELE + 1}},
  };
  child_map child_m = build_child_map(par_map);
  gt_sites sites;
  void SetUp() {
    sites = gt_sites{std::make_shared<LevelGenotypedSite>(),
                     std::make_shared<LevelGenotypedSite>(),
                     std::make_shared<LevelGenotypedSite>()};
  }
};

TEST_F(LevelGenotyperPropagation,
       GivenNestingStructure_CorrectGenotypeNullifying) {
  auto site1 = std::dynamic_pointer_cast<LevelGenotypedSite>(sites.at(1));
  site1->set_num_haplogroups(5);

  // SiteID 9 will get nulled by site 7.
  // Then when site 5 nulls site 7, I want site 9 to signal it is already
  // nulled.
  auto site2 = std::dynamic_pointer_cast<LevelGenotypedSite>(sites.at(2));
  site2->set_num_haplogroups(5);

  LevelGenotyper g{child_m, sites};

  EXPECT_FALSE(site2->is_null());
  // I want site 9 to be invalidated in this call
  g.invalidate_if_needed(7, AlleleIds{1});
  EXPECT_TRUE(site2->is_null());

  EXPECT_FALSE(site1->is_null());
  // And this call to invalidate site 7
  g.invalidate_if_needed(5, AlleleIds{0});
  EXPECT_TRUE(site1->is_null());
}

TEST_F(LevelGenotyperPropagation, CorrectFilterDownPropagation) {
  LevelGenotyper g{child_m, sites};
  g.downpropagate_filter("AMBIG", 5);
  EXPECT_TRUE(sites.at(1)->has_filter("AMBIG"));
  EXPECT_TRUE(sites.at(2)->has_filter("AMBIG"));
}

TEST_F(LevelGenotyperPropagation, CorrectFilterUpPropagation) {
  LevelGenotyper g{child_m, sites};
  sites.at(1)->set_filter("AMBIG");
  g.uppropagate_filter("AMBIG", 5);
  EXPECT_TRUE(sites.at(0)->has_filter("AMBIG"));
}
