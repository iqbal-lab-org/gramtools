/**
 * @file
 * Test high-level quasimapping routine: searching for full kmers or full reads.
 * Assessing results is in terms of SearchStates produced or coverage recorded.
 *
 * Suites:
 *  - SearchStates: test that you produce the right search states
 *  - Coverage: test that mapping increments the right allele sum coverage,
 * grouped allele counts coverage, and/or per base coverage.
 *
 *  A "_Nested" suffix is added for nested PRGs.
 *
 */

#include <stdexcept>

#include "genotype/quasimap/coverage/allele_base.hpp"
#include "genotype/quasimap/quasimap.hpp"
#include "genotype/quasimap/search/BWT_search.hpp"
#include "gtest/gtest.h"
#include "submod_resources.hpp"
#include "test_resources.hpp"

using namespace gram::submods;

TEST(ReverseComplementRead, GivenRead_ReverseComplementReadReturned) {
  gram::Sequence read = {1, 2, 1, 3, 4};
  auto result = gram::reverse_complement_read(read);
  gram::Sequence expected = {1, 2, 4, 3, 4};
  EXPECT_EQ(result, expected);
}

TEST(GetKmers, GivenKmersAtOffsets_CorrectExtractionAndThrowsIfKmerDoesNotFit) {
  auto read = encode_dna_bases("accgaat");
  uint32_t kmer_size = 4;
  std::vector<std::string> expected_fits{"accg", "ccga", "cgaa", "gaat"};
  Sequences encoded_fits;
  for (auto const &kmer : expected_fits)
    encoded_fits.push_back(encode_dna_bases(kmer));
  std::size_t max_offset = read.size() - kmer_size;
  for (std::size_t offset = 0; offset <= max_offset; ++offset) {
    auto result = get_kmer_in_read(kmer_size, offset, read);
    EXPECT_EQ(result, encoded_fits.at(offset));
  }
  EXPECT_THROW(get_kmer_in_read(kmer_size, max_offset + 1, read),
               std::invalid_argument);
}

TEST(GetKmers, GivenReadAndKmerSize_CorrectLastKmerReturned) {
  auto read = encode_dna_bases("accgaatt");
  uint32_t kmer_size = 3;
  auto result = get_last_kmer_in_read(kmer_size, read);
  auto expected = encode_dna_bases("att");
  EXPECT_EQ(result, expected);
}

TEST(KmersAllInRead, GivenKmerIndex_AllKmersInReadMustBeIndexed) {
  uint32_t kmer_size = 4;
  KmerIndex index{{encode_dna_bases("accg"), SearchStates{}},
                  {encode_dna_bases("ccgt"), SearchStates{}}};
  auto read1 = encode_dna_bases("accgt");
  auto read2 = encode_dna_bases("tccgt");
  EXPECT_TRUE(all_read_kmers_occur_in_index(kmer_size, read1, index));
  EXPECT_FALSE(all_read_kmers_occur_in_index(kmer_size, read2, index));
}

TEST(Coverage, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6aG7t8C8CTA");

  const auto read = encode_dna_bases("agccta");

  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 0, 0}, {0, 1}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7t8c8cta");

  const auto read = encode_dna_bases("agtcta");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 0, 0}, {1, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7t8c8cta");

  const auto read = encode_dna_bases("ctgagtcta");

  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 1, 0}, {1, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadCrossTwoSitesAndEndsInSite_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7t8c8cta");

  const auto read = encode_dna_bases("tagtcta");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 0, 1}, {1, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadDoesNotMap_EmptyAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7t8c8cta");

  const auto read = encode_dna_bases("tgtcta");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 0, 0}, {0, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadEndsInAllele_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7t8c8cta");

  const auto read = encode_dna_bases("gctc");

  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{1, 0, 0}, {0, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadStartsInAllele_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6T6AG7T8c8cta");

  const auto read = encode_dna_bases("tagt");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 0, 1}, {1, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadWithNoMatchingKmer_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7t8ta8");

  const auto read = encode_dna_bases("tagc");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 0, 0}, {0, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadMapsToThreePositions_CorrectAlleleCoverage) {
  /**
   * The read has three mapping instances, with two distinct site paths:
   * site 5 only, or site 5 and site 7.
   * Depending on choice of seed, can choose one or the other.
   */
  prg_setup setup;
  setup.setup_numbered_prg("TAG5Tc6g6T6AG7T8c8cta");
  const auto read = encode_dna_bases("tagt");

  // Chooses mapping instance in site 5 only
  SeedSize const random_seed1 = 42;
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats, random_seed1);
  auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{1, 0, 1}, {0, 0}};
  EXPECT_EQ(result, expected);

  // Chooses mapping instance in site 5 + site 7
  SeedSize const random_seed2 = 150;
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats, random_seed2);
  expected = {{1, 0, 2}, {1, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadEntirelyWithinAllele_CoverageRecorded) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5cccc6g6t6ag");

  const auto read = encode_dna_bases("cccc");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{1, 0, 0}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadMapsWithinAllele_SumCoverageIsOne) {
  prg_setup setup;
  setup.setup_numbered_prg("ac5t6cagtagtc6ta");

  Sequence read = encode_dna_bases("gtagt");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 1}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadMapsTwiceWithinAllele_SumCoverageIsOne) {
  prg_setup setup;
  setup.setup_numbered_prg("ac5t6cagtagttttgtagtc6ta");
  setup.parameters.seed = 42;

  Sequence read = encode_dna_bases("gtagt");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 1}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, ReadMapsWithinAlleleAndOutsideSite_CorrectSumCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gtagtac5gtagtact6t6ta");

  SeedSize const random_seed = 29;
  Sequence read = encode_dna_bases("gtagt");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats, random_seed);

  auto const &sumCovResult = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage sumCovExpected = {{1, 0}};
  EXPECT_EQ(sumCovResult, sumCovExpected);

  auto const &pbCovResult =
      coverage::generate::allele_base_non_nested(setup.prg_info);
  SitesAlleleBaseCoverage pbCovExpected{SitePbCoverage{
      PerBaseCoverage{1, 1, 1, 1, 1, 0, 0, 0}, PerBaseCoverage{0}}};
  EXPECT_EQ(pbCovResult, pbCovExpected);
}

TEST(Coverage, ReadEndWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("tac5gta6gtt6ta");

  Sequence read = encode_dna_bases("tacgt");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &sumCovResult = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage sumCovExpected = {{1, 1}};
  EXPECT_EQ(sumCovResult, sumCovExpected);

  auto const &pbCovResult =
      coverage::generate::allele_base_non_nested(setup.prg_info);
  SitesAlleleBaseCoverage pbCovExpected{
      SitePbCoverage{PerBaseCoverage{1, 1, 0}, PerBaseCoverage{1, 1, 0}}};
  EXPECT_EQ(pbCovResult, pbCovExpected);
}

TEST(Coverage, ReadStartWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("c5ccc6agt6ccgt6taa");
  setup.parameters.seed = 39;

  Sequence read = encode_dna_bases("gttaa");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 1, 1}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, EncapsulatedWithinTwoDifferentAlleles_CorrectAlleleSumCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("ac5gtagtact6t6gggtagt6ta");
  setup.parameters.seed = 42;

  Sequence read = encode_dna_bases("gtagt");
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{1, 0, 1}};
  EXPECT_EQ(result, expected);

  auto const &pbCovResult =
      coverage::generate::allele_base_non_nested(setup.prg_info);
  SitesAlleleBaseCoverage pbCovExpected{
      SitePbCoverage{PerBaseCoverage{1, 1, 1, 1, 1, 0, 0, 0},
                     PerBaseCoverage{0}, PerBaseCoverage{0, 0, 1, 1, 1, 1, 1}}};
  EXPECT_EQ(pbCovResult, pbCovExpected);
}

TEST(Coverage, MappingMultipleIdenticalReads_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6T6AG7T8c8cta");

  Sequences reads = {encode_dna_bases("tagt"), encode_dna_bases("tagt")};

  for (const auto &read : reads) {
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                  setup.parameters, setup.quasimap_stats);
  }

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{0, 0, 2}, {2, 0}};
  EXPECT_EQ(result, expected);

  auto const &pbCovResult =
      coverage::generate::allele_base_non_nested(setup.prg_info);
  SitesAlleleBaseCoverage pbCovExpected{
      SitePbCoverage{
          PerBaseCoverage{0},
          PerBaseCoverage{0},
          PerBaseCoverage{2},
      },
      SitePbCoverage{PerBaseCoverage{2}, PerBaseCoverage{0}}};
  EXPECT_EQ(pbCovResult, pbCovExpected);
}

TEST(Coverage, MappingThreeReadsIdenticalKmers_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7t8c8cta");

  Sequences reads = {encode_dna_bases("gagt"), encode_dna_bases("tagt"),
                     encode_dna_bases("cagt")};

  for (const auto &read : reads) {
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                  setup.parameters, setup.quasimap_stats);
  }

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{1, 1, 1}, {3, 0}};
  EXPECT_EQ(result, expected);

  auto const &pbCovResult =
      coverage::generate::allele_base_non_nested(setup.prg_info);
  SitesAlleleBaseCoverage pbCovExpected{
      SitePbCoverage{
          PerBaseCoverage{1},
          PerBaseCoverage{1},
          PerBaseCoverage{1},
      },
      SitePbCoverage{PerBaseCoverage{3}, PerBaseCoverage{0}}};
  EXPECT_EQ(pbCovResult, pbCovExpected);
}

TEST(Coverage, MappingThreeReadsDifferentKmers_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7t8c8cta");

  Sequences reads = {encode_dna_bases("gagt"), encode_dna_bases("tagt"),
                     encode_dna_bases("cagc")};

  for (const auto &read : reads) {
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                  setup.parameters, setup.quasimap_stats);
  }

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{1, 1, 1}, {2, 1}};
  EXPECT_EQ(result, expected);
}

TEST(Coverage, MappingThreeReadsOneReadMapsTwice_CorrectAlleleCoverage) {
  prg_setup setup;
  setup.setup_numbered_prg("gcac5t6g6c6ta7t8c8cta");

  Sequences reads = {
      encode_dna_bases("accta"),
      encode_dna_bases("gcact"),
  };

  SeedSize const random_seed = 200;
  for (const auto &read : reads) {
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                  setup.parameters, setup.quasimap_stats, random_seed);
  }

  const auto &result = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage expected = {{1, 0, 0}, {0, 1}};
  EXPECT_EQ(result, expected);
}

TEST(KmerIndexQuasimap, KmerAbsentFromKmerIndex_NoSearchStatesReturned) {
  auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
  auto prg_info = generate_prg_info(prg_raw);
  Sequence kmer = encode_dna_bases("gtaa");
  Sequences kmers{encode_dna_bases("tagt"), encode_dna_bases("agta"),
                  encode_dna_bases("gtaa")};
  auto kmer_size = 4;
  auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

  auto read = encode_dna_bases("tagtaa");
  auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
  EXPECT_EQ(search_states.size(), 0);
}

TEST(vBWTJump_andBWTExtension, InitiallyInSite_HaveExitedSite) {
  auto prg_raw = encode_prg("gcgct5c6G6t6agtcct");
  auto prg_info = generate_prg_info(prg_raw);
  auto pattern_char = encode_dna_base('t');

  SearchState initial_search_state = {
      SA_Interval{10, 10},  // Starting at char 'g' at index 8 in prg
      VariantSitePath{},
      VariantSitePath{},
  };
  SearchStates initial_search_states = {initial_search_state};

  auto final_search_states = process_read_char_search_states(
      pattern_char, initial_search_states, prg_info);

  EXPECT_EQ(final_search_states.size(), 1);
  const auto &result = final_search_states.front().traversed_path;
  VariantSitePath expected = {
      VariantLocus{5, FIRST_ALLELE + 1},
  };
  EXPECT_EQ(result, expected);
}

class SearchStates_and_Coverage_EndInSite : public ::testing::Test {
 protected:
  void SetUp() { setup.setup_numbered_prg("gcgct5c6g6T6AGTCCt"); }
  Sequence kmer = encode_dna_bases("cc");
  prg_setup setup;
  Sequence read = encode_dna_bases("tagtcc");
};

TEST_F(SearchStates_and_Coverage_EndInSite, MapOneRead_CorrectSearchState) {
  auto search_states =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
  EXPECT_EQ(search_states.size(), 1);

  // Do we end up in right place in SA index?
  auto search_state = search_states.front();
  auto result = search_state.sa_interval;
  SA_Interval expected = {14, 14};
  EXPECT_EQ(result, expected);

  auto path_result = search_state.traversing_path;
  VariantSitePath path_expected = {VariantLocus{5, ALLELE_UNKNOWN}};
  EXPECT_EQ(path_result, path_expected);
}

TEST_F(SearchStates_and_Coverage_EndInSite, MapOneRead_CorrectCoverage) {
  quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &AlSumResult = setup.coverage.allele_sum_coverage;
  AlleleSumCoverage AlSumExpected = {{0, 0, 1}};
  EXPECT_EQ(AlSumResult, AlSumExpected);

  auto const &pbCovResult =
      coverage::generate::allele_base_non_nested(setup.prg_info);
  SitesAlleleBaseCoverage pbCovExpected{SitePbCoverage{
      PerBaseCoverage{0}, PerBaseCoverage{0}, PerBaseCoverage{1}}};
  EXPECT_EQ(pbCovResult, pbCovExpected);
}

TEST(SearchStates, StartInSiteAndMapOut_CorrectVarLocusPath) {
  prg_setup setup;
  setup.setup_numbered_prg("gcGCT5C6g6t6agtcct");

  auto read = encode_dna_bases("gcgctc");
  Sequence kmer = encode_dna_bases("tc");
  auto search_states =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
  EXPECT_EQ(search_states.size(), 1);

  auto result = search_states.front().traversed_path;
  VariantSitePath expected = {VariantLocus{5, FIRST_ALLELE}};
  EXPECT_EQ(result, expected);
}

TEST(SearchStates, StartOutOfSiteAndMapThrough_CorrectVarLocusPath) {
  prg_setup setup;
  setup.setup_numbered_prg("gcgcT5c6G6t6AGtcct");

  auto read = encode_dna_bases("gctgag");
  Sequence kmer = encode_dna_bases("ag");
  auto search_states =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);

  EXPECT_EQ(search_states.size(), 1);

  auto result = search_states.front().traversed_path;
  VariantSitePath expected = {VariantLocus{5, FIRST_ALLELE + 1}};
  EXPECT_EQ(result, expected);
}

TEST(SearchStates, ReadCrossingTwoAlleles_CorrectVarLocusPaths) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7T8c8CT");

  auto read = encode_dna_bases("cagtct");
  Sequence kmer = encode_dna_bases("ct");
  auto search_states =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
  EXPECT_EQ(search_states.size(), 1);

  auto traversed_path = search_states.front().traversed_path;
  VariantSitePath expected_traversed = {VariantLocus{7, FIRST_ALLELE}};
  EXPECT_EQ(traversed_path, expected_traversed);

  auto traversing_path = search_states.front().traversing_path;
  VariantSitePath expected_traversing = {VariantLocus{5, ALLELE_UNKNOWN}};
  EXPECT_EQ(traversing_path, expected_traversing);
}

TEST(SearchStates, StartWithinAlleleEndWithinAnother_CorrectVarLocusPath) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5c6g6t6ag7GAG8c8ct");

  auto read = encode_dna_bases("caggag");
  Sequence kmer = encode_dna_bases("ag");
  auto search_states =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
  EXPECT_EQ(search_states.size(), 1);

  auto traversed_path = search_states.front().traversed_path;
  VariantSitePath expected_traversed = {VariantLocus{7, FIRST_ALLELE}};
  EXPECT_EQ(traversed_path, expected_traversed);

  auto traversing_path = search_states.front().traversing_path;
  VariantSitePath expected_traversing = {VariantLocus{5, ALLELE_UNKNOWN}};
  EXPECT_EQ(traversing_path, expected_traversing);
}

/*
 * A case where we end the read mapping inside several alleles of the same site.
 * We test: correct indexing, correct base extension, correct allele id
 * specification.
 */
TEST(MultiStep, RunIndexingExtensionIdSpecification_CorrectOutputs) {
  prg_setup setup;
  setup.setup_numbered_prg("gct5gC6aC6C6t6Cg", 1);

  // We expect five occurrences of 'C' at this stage, in a single SA interval
  Sequence kmer = encode_dna_bases("c");
  auto search_states = setup.kmer_index.at(kmer);
  EXPECT_EQ(search_states.size(), 1);
  SA_Interval sa = search_states.front().sa_interval;
  EXPECT_EQ(sa.second - sa.first + 1, 5);

  // Next up, look for a C
  int_Base pattern_char = 2;
  search_states = process_read_char_search_states(pattern_char, search_states,
                                                  setup.prg_info);

  // concurrent allele querying
  // Expect three occurrences of 'CC' at this stage, in a single SA interval
  EXPECT_EQ(search_states.size(), 1);
  EXPECT_EQ(search_states.front().traversing_path.back().second,
            ALLELE_UNKNOWN);
  EXPECT_EQ(search_states.front().sa_interval.second -
                search_states.front().sa_interval.first + 1,
            3);
}

TEST(SearchStates, OneMappingEncapsulatedByAllele) {
  prg_setup setup;
  setup.setup_numbered_prg("t5c6gCTTAGT6aa");

  auto read = encode_dna_bases("cttagt");
  Sequence kmer = encode_dna_bases("gt");
  auto search_states =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
  EXPECT_EQ(search_states.size(), 1);

  VariantLocus cov = {5, FIRST_ALLELE + 1};
  EXPECT_EQ(search_states.front().traversed_path.front(), cov);
}

TEST(SearchStates, StartAndEndInSite_CorrectSearchStates) {
  prg_setup setup;
  setup.setup_numbered_prg("t5c6gcttagtacgcttagt6aa");

  auto read = encode_dna_bases("cttagt");
  Sequence kmer = encode_dna_bases("gt");
  auto result =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);

  SearchStates expected = {SearchState{
      SA_Interval{7, 8},
      VariantSitePath{VariantLocus{5, FIRST_ALLELE + 1}},
      VariantSitePath{},
  }};

  EXPECT_EQ(result, expected);
}

TEST(SearchStates_Nested, MapIntoAndOutOfNestedSite_CorrectSearchStates) {
  prg_setup setup;
  setup.setup_bracketed_prg("a[c,g[ct,t]a]c");

  auto read = encode_dna_bases("agtac");
  Sequence kmer = encode_dna_bases("ac");
  auto result =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);

  SearchStates expected = {SearchState{
      SA_Interval{1, 1},
      VariantSitePath{VariantLocus{7, FIRST_ALLELE + 1},
                      VariantLocus{5, FIRST_ALLELE + 1}},
      VariantSitePath{},
  }};
  EXPECT_EQ(result, expected);
}

/*
PRG: T[A[C,G][C,G],]T
i	BWT	SA	text_suffix
0	T	16	0
1	5	2	A 7 C 8 G 8 9 C 10 G 10 6 6 T 0
2	7	4	C 8 G 8 9 C 10 G 10 6 6 T 0
3	9	9	C 10 G 10 6 6 T 0
4	8	6	G 8 9 C 10 G 10 6 6 T 0
5	10	11	G 10 6 6 T 0
6	6	15	T 0
7	0	0	T 5 A 7 C 8 G 8 9 C 10 G 10 6 6 T 0
8	T	1	5 A 7 C 8 G 8 9 C 10 G 10 6 6 T 0
9	6	14	6 T 0
10	10	13	6 6 T 0
11	A	3	7 C 8 G 8 9 C 10 G 10 6 6 T 0
12	C	5	8 G 8 9 C 10 G 10 6 6 T 0
13	G	7	8 9 C 10 G 10 6 6 T 0
14	8	8	9 C 10 G 10 6 6 T 0
15	C	10	10 G 10 6 6 T 0
16	G	12	10 6 6 T 0
*/
TEST(ReadQuasimap_Nested, MapThroughDeletionAndExitEntry_CorrectSearchStates) {
  prg_setup setup;
  setup.setup_bracketed_prg("t[a[c,g][c,g],]t", 1);

  auto read = encode_dna_bases("tt");
  Sequence kmer = encode_dna_bases("t");
  auto result_direct_deletion =
      search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);

  SearchStates expected_direct_deletion = {SearchState{
      SA_Interval{7, 7},
      VariantSitePath{VariantLocus{5, FIRST_ALLELE + 1}},
      VariantSitePath{},
  }};
  EXPECT_EQ(result_direct_deletion, expected_direct_deletion);

  auto read2 = encode_dna_bases("tacct");
  auto result_exit_entry =
      search_read_backwards(read2, kmer, setup.kmer_index, setup.prg_info);

  SearchStates expected_exit_entry = {SearchState{
      SA_Interval{7, 7},
      VariantSitePath{VariantLocus{9, FIRST_ALLELE},
                      VariantLocus{7, FIRST_ALLELE},
                      VariantLocus{5, FIRST_ALLELE}},
      VariantSitePath{},
  }};
  EXPECT_EQ(result_exit_entry, expected_exit_entry);
}

class Coverage_Nested_DoubleNesting : public ::testing::Test {
  // Double nesting meaning a bubble inside a bubble inside a bubble
 protected:
  void SetUp() { setup.setup_bracketed_prg("A[[A[CCC,c],t],g]TA"); }
  prg_setup setup;
  prg_positions positions{
      0, 3, 5, 9, 12, 15, 17};  // All the nodes in the cov graph with sequence

  Sequence read1 = encode_dna_bases("AACCCTA");
  Sequence read2 = encode_dna_bases("CTA");
};

TEST_F(Coverage_Nested_DoubleNesting,
       ReadEndsInsideNestedSite_CorrectCoverage) {
  // PRG: "A[[A[CCC,c],t],g]TA"; Read: "AACCCTA"
  quasimap_read(read1, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
  // The read is compatible with the first allele of all three sites in the PRG
  SitesGroupedAlleleCounts expectedGpAlCounts = {
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
  };
  EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

  auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
  SitePbCoverage expectedPbCov{PerBaseCoverage{},        PerBaseCoverage{1},
                               PerBaseCoverage{1, 1, 1}, PerBaseCoverage{0},
                               PerBaseCoverage{0},       PerBaseCoverage{0},
                               PerBaseCoverage{}};
  EXPECT_EQ(PbCov, expectedPbCov);
}

TEST_F(Coverage_Nested_DoubleNesting, ReadMultiMaps_CorrectCoverage) {
  // PRG: "A[[A[CCC,c],t],g]TA"; Read: "CTA"
  quasimap_read(read2, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
  // The read is compatible with the two alleles of the most nested site in the
  // PRG string
  SitesGroupedAlleleCounts expectedGpAlCounts = {
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
      GroupedAlleleCounts{{AlleleIds{0, 1}, 1}},
  };
  EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

  auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
  SitePbCoverage expectedPbCov{PerBaseCoverage{},        PerBaseCoverage{0},
                               PerBaseCoverage{0, 0, 1}, PerBaseCoverage{1},
                               PerBaseCoverage{0},       PerBaseCoverage{0},
                               PerBaseCoverage{}};
  EXPECT_EQ(PbCov, expectedPbCov);
}

class Coverage_Nested_SingleNestingPlusSNP : public ::testing::Test {
 protected:
  void SetUp() { setup.setup_bracketed_prg("a[t[tt,t]t,a[at,]a]g[c,g]"); }
  prg_setup setup;
  prg_positions positions{
      0,  2,  4,  7,  9, 11,
      13, 17, 19, 21, 23};  // All the nodes in the cov graph with sequence

  Sequence read1 = encode_dna_bases("ATTTTGC");
  Sequence read2 = encode_dna_bases("TT");
  Sequence read3 = encode_dna_bases("AAAGG");
};

TEST_F(Coverage_Nested_SingleNestingPlusSNP,
       FullyCrossingRead_CorrectCoverage) {
  // PRG: "A[T[TT,T]T,a[at,]a]G[C,g]" ; Read: "ATTTTGC"
  quasimap_read(read1, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
  // The read is compatible with the two alleles of the most nested site in the
  // PRG string
  SitesGroupedAlleleCounts expectedGpAlCounts = {
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
      GroupedAlleleCounts{},
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
  };
  EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

  auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
  SitePbCoverage expectedPbCov{
      PerBaseCoverage{},     PerBaseCoverage{1}, PerBaseCoverage{1, 1},
      PerBaseCoverage{0},    PerBaseCoverage{1}, PerBaseCoverage{0},
      PerBaseCoverage{0, 0}, PerBaseCoverage{0}, PerBaseCoverage{},
      PerBaseCoverage{1},    PerBaseCoverage{0}};
  EXPECT_EQ(PbCov, expectedPbCov);
}

TEST_F(Coverage_Nested_SingleNestingPlusSNP,
       VeryMultiMappingRead_CorrectCoverage) {
  // PRG: "A[T[TT,T]T,a[at,]a]G[C,g]" ; Read: "TT"
  // This read should have 5 mapping instances: one is encapsulated(=empty
  // traversing and traversed), two are in 'traversing' mode, two are in
  // 'traversed' mode. All are encapsulated inside site 0 as well.

  quasimap_read(read2, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
  // The read is compatible with the two alleles of the most nested site in the
  // PRG string
  SitesGroupedAlleleCounts expectedGpAlCounts = {
      GroupedAlleleCounts{{AlleleIds{0}, 1}},
      GroupedAlleleCounts{{AlleleIds{0, 1}, 1}},
      GroupedAlleleCounts{},
      GroupedAlleleCounts{},
  };
  EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

  auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
  SitePbCoverage expectedPbCov{
      PerBaseCoverage{},     PerBaseCoverage{1}, PerBaseCoverage{1, 1},
      PerBaseCoverage{1},    PerBaseCoverage{1}, PerBaseCoverage{0},
      PerBaseCoverage{0, 0}, PerBaseCoverage{0}, PerBaseCoverage{},
      PerBaseCoverage{0},    PerBaseCoverage{0}};
  EXPECT_EQ(PbCov, expectedPbCov);
}

TEST_F(Coverage_Nested_SingleNestingPlusSNP,
       MapThroughDirectDeletion_CorrectCoverage) {
  // PRG: "A[t[tt,t]t,A[at,]A]G[c,G]" ; Read: "AAAGG"
  quasimap_read(read3, setup.coverage, setup.kmer_index, setup.prg_info,
                setup.parameters, setup.quasimap_stats);

  const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
  SitesGroupedAlleleCounts expectedGpAlCounts = {
      GroupedAlleleCounts{{AlleleIds{1}, 1}},
      GroupedAlleleCounts{},
      GroupedAlleleCounts{{AlleleIds{1}, 1}},
      GroupedAlleleCounts{{AlleleIds{1}, 1}},
  };
  EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

  auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
  SitePbCoverage expectedPbCov{
      PerBaseCoverage{},     PerBaseCoverage{0}, PerBaseCoverage{0, 0},
      PerBaseCoverage{0},    PerBaseCoverage{0}, PerBaseCoverage{1},
      PerBaseCoverage{0, 0}, PerBaseCoverage{1}, PerBaseCoverage{},
      PerBaseCoverage{0},    PerBaseCoverage{1}};
  EXPECT_EQ(PbCov, expectedPbCov);
}
