#include "genotype/read_stats.hpp"
#include "gtest/gtest.h"
#include "mocks.hpp"
#include "prg/coverage_graph.hpp"
#include "submod_resources.hpp"
#include "test_resources.hpp"

using namespace gram;

/**
 * Per-base error rate
 */

TEST(ReadProcessingStats, GivenTwoGenomicReads_CorrectStats) {
  GenomicRead_vector reads{
      GenomicRead{"Read1", "AAAA", "5555"},  // 5: ASCII 53, Q-score 20 (Phred
                                             // +33 scale), error prob 0.01.
      GenomicRead{"Read1", "TTTT", "5555"}};

  ReadStats r;
  r.compute_base_error_rate(reads);

  EXPECT_EQ(r.get_num_bases_processed(), 8);
  EXPECT_EQ(r.get_max_read_len(), 4);
  EXPECT_FLOAT_EQ(r.get_mean_pb_error(), 0.01);
}

TEST(ReadProcessingStats, GivenOneOKAndOneEmptyGenomicRead_CorrectStats) {
  GenomicRead_vector reads{
      GenomicRead{"Read1", "AAA",
                  "???"},  // ASCII 63, Q-score 30, error prob 0.001
      GenomicRead{"Read1", "", ""}};

  ReadStats r;
  r.compute_base_error_rate(reads);

  EXPECT_EQ(r.get_num_no_qual_reads(), 1);
  EXPECT_FLOAT_EQ(r.get_mean_pb_error(), 0.001);
}

/**
 * Coverage mean and variance
 * Notes:
 *   - coverage statistics only measured from level 1 sites (sites
 *     not nested in any others) to avoid double-counting
 *   - for each site, allele with max coverage is extracted,
 *     and per-base coverage computed from it
 */

TEST(MaxHaplogroupCoverage, GivenGroupedAlleleCoverage_CorrectMax) {
  using res_cov = ReadStats::haplogroup_cov;
  auto result = ReadStats::get_max_cov_haplogroup({});
  auto expected = res_cov{0, 0};
  EXPECT_EQ(result, expected);

  GroupedAlleleCounts gped_covs{
      {AlleleIds{0, 1}, 2}, {AlleleIds{0}, 3}, {AlleleIds{1}, 4}};
  // The single allele with max coverage on it is returned along with
  // the total coverage on the site
  result = ReadStats::get_max_cov_haplogroup(gped_covs);
  expected = res_cov{1, 9};
}

using namespace gram::submods;
class TestReadMappingStats : public ::testing::Test {
 protected:
  void SetUp() {
    auto prg = prg_string_to_ints("[AC[T,G]AC,GT[A,T]T]A[AA,C]T");
    cov_graph = coverage_Graph(prg);
  }

  Coverage cov{
      {},
      SitesGroupedAlleleCounts{
          GroupedAlleleCounts{{AlleleIds{1}, 60}},
          GroupedAlleleCounts{{AlleleIds{1}, 2}, {AlleleIds{0}, 1}},
          GroupedAlleleCounts{{AlleleIds{0}, 19}, {AlleleIds{0, 1}, 1}},
          GroupedAlleleCounts{}},
      {}};

  coverage_Graph cov_graph;
  ReadStats stats;
};

TEST_F(TestReadMappingStats, ExtractMaxCovAlleleSite1) {
  auto site_pair = get_bubble_nodes(cov_graph.bubble_map, 7);
  auto extracted = stats.extract_max_coverage_allele(
      cov.grouped_allele_counts, site_pair.first, site_pair.second);
  EXPECT_EQ(extracted.first.sequence, "G");
  EXPECT_EQ(extracted.second, 2);
}

TEST_F(TestReadMappingStats, ExtractMaxCovAlleleSite2) {
  auto site_pair = get_bubble_nodes(cov_graph.bubble_map, 9);
  auto extracted = stats.extract_max_coverage_allele(
      cov.grouped_allele_counts, site_pair.first, site_pair.second);
  EXPECT_EQ(extracted.first.sequence, "A");
  EXPECT_EQ(extracted.second, 20);
}

TEST_F(TestReadMappingStats, ExtractMaxCovAlleleSite3) {
  auto site_pair = get_bubble_nodes(cov_graph.bubble_map, 11);
  auto extracted = stats.extract_max_coverage_allele(
      cov.grouped_allele_counts, site_pair.first, site_pair.second);
  EXPECT_EQ(extracted.first.sequence, "AA");
  EXPECT_EQ(extracted.second, 0);
}

TEST_F(TestReadMappingStats, ExtractMaxCovAlleleSite0) {
  auto site_pair = get_bubble_nodes(cov_graph.bubble_map, 5);
  auto extracted = stats.extract_max_coverage_allele(
      cov.grouped_allele_counts, site_pair.first, site_pair.second);
  EXPECT_EQ(extracted.first.sequence, "GTAT");
  EXPECT_EQ(extracted.second, 60);
}

TEST(TestMeanAndVarCovComputation, GivenMockReturnedAlleles_CorrectStats) {
  using ::testing::_;
  using ::testing::Return;
  using returned = AbstractReadStats::allele_and_cov;
  MockReadStats stats;
  auto prg = prg_string_to_ints("A[A,T]C[T,C]");
  coverage_Graph cov_graph{prg};  // Only used for identifying var sites

  EXPECT_CALL(stats, extract_max_coverage_allele(_, _, _))
      .Times(2)
      .WillOnce(Return(returned{Allele{"AT", {10, 20}}, 20}))
      .WillOnce(Return(returned{Allele{"", {}}, 5}));
  stats.compute_coverage_depth(Coverage{}, cov_graph);

  // Expect the mean of 15 ((10 + 20) / 2) and 5 (direct deletion)
  EXPECT_DOUBLE_EQ(stats.get_mean_cov(), 10);
  EXPECT_DOUBLE_EQ(stats.get_var_cov(), 25);
}

/**
 * Integration tests. Rely on proper mapping, coverage recording
 * and parental map formation so also test those.
 */
TEST(ReadMappingStats,
     GivenThreeMappedReadsNonNestedPRG_CorrectMappingRelatedStats) {
  GenomicRead_vector reads{
      GenomicRead{"Read1", "AAA", "###"},  // '#' = Q-score of 2
      GenomicRead{"Read2", "AAA", "###"},
      GenomicRead{"Read3", "GCAAA", "#####"},
      GenomicRead{"Read4", "GCAAA", "#####"}};

  prg_setup setup;
  Sequences kmers{encode_dna_bases("AA")};
  setup.setup_numbered_prg("G5CAAA6AA6T7G8C8GGG", kmers);
  setup.quasimap_reads(reads);

  auto stats = setup.read_stats;
  // Map 4 reads to site 1, and 0 to site 2.
  // Estimated per base cov at site 1 is (2 + 4 + 4 + 4) / 4 = 3.5
  // (allele 'AAAA' is single most supported). Mean is (3.5 + 0) / 2
  EXPECT_DOUBLE_EQ(stats.get_mean_cov(), 1.75);
  EXPECT_DOUBLE_EQ(stats.get_var_cov(), 3.0625);
  EXPECT_EQ(stats.get_num_sites_noCov(), 1);
  EXPECT_EQ(stats.get_num_sites_total(), 2);
}

TEST(ReadMappingStats,
     GivenTwoMappedReadsNestedPRG_CorrectMappingRelatedStats) {
  GenomicRead_vector reads{
      GenomicRead{"Read1", "GGGGGCCC", "IIIIIIIII"},  // 'I' = Q-score of 40
      GenomicRead{"Read2", "GCCCC", "IIII"},
      GenomicRead{"Read3", "GCCCC", "IIII"},
      GenomicRead{"Read4", "GCCC", "IIII"},
  };  // Read4 compatible with both alleles of parent site

  prg_setup setup;
  Sequences kmers{encode_dna_bases("CC")};
  setup.setup_bracketed_prg("G[GG[G,A]G,C]CCC", kmers);
  setup.quasimap_reads(reads);

  auto stats = setup.read_stats;
  // 3: because single allele with most cov is 'C' in parent site, and 3 reads
  // go through it
  EXPECT_DOUBLE_EQ(stats.get_mean_cov(), 3);
  EXPECT_DOUBLE_EQ(stats.get_var_cov(), 0);
  EXPECT_EQ(stats.get_num_sites_noCov(), 0);
  EXPECT_EQ(stats.get_num_sites_total(), 1);
}
