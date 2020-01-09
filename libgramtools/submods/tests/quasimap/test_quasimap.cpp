/**
 * @file
 * Test high-level quasimapping routine: searching for full kmers or full reads.
 * Assessing results is in terms of SearchStates produced or coverage recorded.
 *
 * Suites:
 *  - SearchStates: test that you produce the right search states
 *  - Coverage: test that mapping increments the right allele sum coverage, grouped allele counts coverage, and/or
 *      per base coverage.
 *
 *  A "_Nested" suffix is added for nested PRGs.
 *
 */

#include "gtest/gtest.h"

#include "src_common/generate_prg.hpp"
#include "kmer_index/build.hpp"
#include "genotype/quasimap/coverage/common.hpp"
#include "genotype/quasimap/quasimap.hpp"
#include "common/utils.hpp"
#include "genotype/quasimap/search/BWT_search.hpp"
#include "genotype/quasimap/coverage/allele_base.hpp"
#include "tests/common.hpp"

using namespace gram;

class prg_setup{
public:
    PRG_Info prg_info;
    Coverage coverage;
    Parameters parameters;
    KmerIndex kmer_index;

    explicit prg_setup() {};
    void setup(std::string raw_prg,
               Sequences kmers){
       auto encoded_prg = encode_prg(raw_prg);
       internal_setup(encoded_prg, kmers);
    }

    void setup_nested(std::string raw_prg,
                      Sequences kmers){
        auto encoded_prg = prg_string_to_ints(raw_prg);
        internal_setup(encoded_prg, kmers);
    }

private:
    void internal_setup(marker_vec encoded_prg, Sequences kmers){
        size_t kmer_size = kmers.front().size();
        for (auto const& kmer : kmers) assert(kmer_size == kmer.size());

        // TODO: the calls to rank_support setup in `generate_prg_info` get somehow lost when leaving its scope
        // and we need to call `init_support`, or rank_support again, in this scope for it to work
        prg_info = generate_prg_info(encoded_prg);

        sdsl::util::init_support(prg_info.rank_bwt_a, &prg_info.dna_bwt_masks.mask_a);
        sdsl::util::init_support(prg_info.rank_bwt_c, &prg_info.dna_bwt_masks.mask_c);
        sdsl::util::init_support(prg_info.rank_bwt_g, &prg_info.dna_bwt_masks.mask_g);
        sdsl::util::init_support(prg_info.rank_bwt_t, &prg_info.dna_bwt_masks.mask_t);

        sdsl::util::init_support(prg_info.prg_markers_rank, &prg_info.prg_markers_mask);
        sdsl::util::init_support(prg_info.prg_markers_select, &prg_info.prg_markers_mask);

        coverage = coverage::generate::empty_structure(prg_info);

        parameters.kmers_size = kmer_size;
        kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);
    }
};

TEST(ReverseComplementRead, GivenRead_ReverseComplementReadReturned) {
    gram::Sequence read = {1, 2, 1, 3, 4};
    auto result = gram::reverse_complement_read(read);
    gram::Sequence expected = {1, 2, 4, 3, 4};
    EXPECT_EQ(result, expected);
}


TEST(GetKmer, GivenReadAndKmerSize_CorrectKmerReturned) {
    auto read = encode_dna_bases("accgaatt");
    uint32_t kmer_size = 3;
    auto result = get_kmer_from_read(kmer_size, read);
    auto expected = encode_dna_bases("att");
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("gccta");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6aG7t8C8CTA", kmers);

    const auto read = encode_dna_bases("agccta");

    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("gtcta");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta",kmers);

    const auto read = encode_dna_bases("agtcta");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("gtcta");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta",kmers);

    const auto read = encode_dna_bases("ctgagtcta");

    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 0},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadCrossTwoSitesAndEndsInSite_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("gtcta");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta",kmers);

    const auto read = encode_dna_bases("tagtcta");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 1},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadDoesNotMap_EmptyAlleleCoverage) {
    Sequence kmer = encode_dna_bases("gtcta");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta",kmers);

    const auto read = encode_dna_bases("tgtcta");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadEndsInAllele_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("ctc");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta", kmers);

    const auto read = encode_dna_bases("gctc");

    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadStartsInAllele_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("agt");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6T6AG7T8c8cta",kmers);

    const auto read = encode_dna_bases("tagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 1},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadWithNoMatchingKmer_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("agt");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta", kmers);

    const auto read = encode_dna_bases("tagc");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadMapsToThreePositions_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("agt");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("tag5tc6g6t6ag7t8c8cta",kmers);

    setup.parameters.seed = 42;
    const auto read = encode_dna_bases("tagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadEntirelyWithinAllele_CoverageRecorded) {
    Sequence kmer = encode_dna_bases("ccc");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5cccc6g6t6ag",kmers);

    const auto read = encode_dna_bases("cccc");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadMapsWithinAllele_SumCoverageIsOne) {
    Sequences kmers = {
            encode_dna_bases("agt"),
    };
    prg_setup setup;
    setup.setup("ac5t6cagtagtc6ta",kmers);

    Sequence read = encode_dna_bases("gtagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadMapsTwiceWithinAllele_SumCoverageIsOne) {
    Sequences kmers = {
            encode_dna_bases("agt"),
    };
    prg_setup setup;
    setup.setup("ac5t6cagtagttttgtagtc6ta",kmers);
    setup.parameters.seed = 42;

    Sequence read = encode_dna_bases("gtagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, ReadMapsWithinAlleleAndOutsideSite_CorrectSumCoverage) {
    Sequences kmers = {
            encode_dna_bases("agt"),
    };
    prg_setup setup;
    setup.setup("gtagtac5gtagtact6t6ta",kmers);
    setup.parameters.seed = 39;

    Sequence read = encode_dna_bases("gtagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    auto const& sumCovResult = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage sumCovExpected = {
            {1, 0}
    };
    EXPECT_EQ(sumCovResult, sumCovExpected);

    auto const& pbCovResult = coverage::generate::allele_base_non_nested(setup.prg_info);
    SitesAlleleBaseCoverage pbCovExpected{
            AlleleCoverage{
                    BaseCoverage{1, 1, 1, 1, 1, 0, 0, 0}, BaseCoverage{0}
            }
    };
    EXPECT_EQ(pbCovResult, pbCovExpected);
}


TEST(Coverage, ReadEndWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
    Sequences kmers = {
            encode_dna_bases("cgt"),
    };
    prg_setup setup;
    setup.setup("tac5gta6gtt6ta", kmers);

    Sequence read = encode_dna_bases("tacgt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &sumCovResult = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage sumCovExpected = {
            {1, 1}
    };
    EXPECT_EQ(sumCovResult, sumCovExpected);

    auto const& pbCovResult = coverage::generate::allele_base_non_nested(setup.prg_info);
    SitesAlleleBaseCoverage pbCovExpected{
       AlleleCoverage{
           BaseCoverage{1, 1, 0}, BaseCoverage{1, 1, 0}
       }
    };
    EXPECT_EQ(pbCovResult, pbCovExpected);
}


TEST(Coverage, ReadStartWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
    Sequences kmers = {
            encode_dna_bases("taa"),
    };
    prg_setup setup;
    setup.setup("c5ccc6agt6ccgt6taa", kmers);
    setup.parameters.seed = 39;

    Sequence read = encode_dna_bases("gttaa");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, EncapsulatedWithinTwoDifferentAlleles_CorrectAlleleSumCoverage) {
    Sequences kmers = {
            encode_dna_bases("agt"),
    };
    prg_setup setup;
    setup.setup("ac5gtagtact6t6gggtagt6ta", kmers);
    setup.parameters.seed = 42;

    Sequence read = encode_dna_bases("gtagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1}
    };
    EXPECT_EQ(result, expected);

    auto const& pbCovResult = coverage::generate::allele_base_non_nested(setup.prg_info);
    SitesAlleleBaseCoverage pbCovExpected{
            AlleleCoverage{
                    BaseCoverage{1, 1, 1, 1, 1, 0, 0, 0}, BaseCoverage{0},
                    BaseCoverage{0, 0, 1, 1, 1, 1, 1}
            }
    };
    EXPECT_EQ(pbCovResult, pbCovExpected);
}


TEST(Coverage, MappingMultipleIdenticalReads_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("agt");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6T6AG7T8c8cta", kmers);

    Sequences reads = {
            encode_dna_bases("tagt"),
            encode_dna_bases("tagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);
    }

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 2},
            {2, 0}
    };
    EXPECT_EQ(result, expected);

    auto const& pbCovResult = coverage::generate::allele_base_non_nested(setup.prg_info);
    SitesAlleleBaseCoverage pbCovExpected{
            AlleleCoverage{
                    BaseCoverage{0}, BaseCoverage{0}, BaseCoverage{2},
            },
            AlleleCoverage{
                    BaseCoverage{2}, BaseCoverage{0}
            }
    };
    EXPECT_EQ(pbCovResult, pbCovExpected);
}

TEST(Coverage, MappingThreeReadsIdenticalKmers_CorrectAlleleCoverage) {
    Sequence kmer = encode_dna_bases("agt");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta", kmers);

    Sequences reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt"),
            encode_dna_bases("cagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);
    }

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 1, 1},
            {3, 0}
    };
    EXPECT_EQ(result, expected);

    auto const& pbCovResult = coverage::generate::allele_base_non_nested(setup.prg_info);
    SitesAlleleBaseCoverage pbCovExpected{
            AlleleCoverage{
                    BaseCoverage{1}, BaseCoverage{1}, BaseCoverage{1},
            },
            AlleleCoverage{
                    BaseCoverage{3}, BaseCoverage{0}
            }
    };
    EXPECT_EQ(pbCovResult, pbCovExpected);
}


TEST(Coverage, MappingThreeReadsDifferentKmers_CorrectAlleleCoverage) {
    Sequences kmers = {
            encode_dna_bases("agt"),
            encode_dna_bases("agc"),
    };
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta", kmers);

    Sequences reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt"),
            encode_dna_bases("cagc")
    };

    for (const auto &read: reads) {
        quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);
    }

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 1, 1},
            {2, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Coverage, MappingThreeReadsOneReadMapsTwice_CorrectAlleleCoverage) {
    Sequences kmers = {
            encode_dna_bases("cta"),
            encode_dna_bases("act"),
    };
    prg_setup setup;
    setup.setup("gcac5t6g6c6ta7t8c8cta", kmers);
    setup.parameters.seed = 42;

    Sequences reads = {
            encode_dna_bases("accta"),
            encode_dna_bases("gcact"),
    };

    for (const auto &read: reads) {
        quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);
    }

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(KmerIndexQuasimap, KmerAbsentFromKmerIndex_NoSearchStatesReturned) {
    auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
    auto prg_info = generate_prg_info(prg_raw);

    Sequence kmer = encode_dna_bases("gtaa");
    Sequences kmers = {kmer};
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
            SA_Interval{10, 10}, //Starting at char 'g' at index 8 in prg
            VariantSitePath{},
            VariantSitePath{},
            SearchVariantSiteState::unknown,
    };
    SearchStates initial_search_states = {initial_search_state};

    auto final_search_states = process_read_char_search_states(pattern_char, initial_search_states, prg_info);

    EXPECT_EQ(final_search_states.size(), 1);
    const auto &result = final_search_states.front().traversed_path;
    VariantSitePath expected = {
            VariantLocus{5, 2},
    };
    EXPECT_EQ(result, expected);
}

class SearchStates_and_Coverage_EndInSite : public::testing::Test{
protected:
    void SetUp(){
        kmer = encode_dna_bases("gtcc");
        Sequences kmers = {kmer};
        setup.setup("gcgct5c6g6T6AGTCCt", kmers);
    }
    Sequence kmer;
    prg_setup setup;
    Sequence read = encode_dna_bases("tagtcc");
};

TEST_F(SearchStates_and_Coverage_EndInSite, MapOneRead_CorrectSearchState) {
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
    EXPECT_EQ(search_states.size(), 1);

    //Do we end up in right place in SA index?
    auto search_state = search_states.front();
    auto result = search_state.sa_interval;
    SA_Interval expected = {14, 14};
    EXPECT_EQ(result, expected);

    auto path_result = search_state.traversing_path;
    VariantSitePath path_expected = {
            VariantLocus{5, ALLELE_UNKNOWN}
    };
    EXPECT_EQ(path_result, path_expected);
}

TEST_F(SearchStates_and_Coverage_EndInSite, MapOneRead_CorrectCoverage) {
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &AlSumResult = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage AlSumExpected = {
            {0, 0, 1}
    };
    EXPECT_EQ(AlSumResult, AlSumExpected );

    auto const& pbCovResult = coverage::generate::allele_base_non_nested(setup.prg_info);
    SitesAlleleBaseCoverage pbCovExpected{
            AlleleCoverage{
                    BaseCoverage{0}, BaseCoverage{0}, BaseCoverage{1}
            }
    };
    EXPECT_EQ(pbCovResult, pbCovExpected);
}

TEST(SearchStates, StartInSiteAndMapOut_CorrectVarLocusPath) {
    Sequence kmer = encode_dna_bases("gctc");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gcGCT5C6g6t6agtcct", kmers);

    auto read = encode_dna_bases("gcgctc");
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto result = search_states.front().traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(SearchStates, StartOutOfSiteAndMapThrough_CorrectVarLocusPath) {
    Sequence kmer = encode_dna_bases("tgag");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gcgcT5c6G6t6AGtcct", kmers);

    auto read = encode_dna_bases("gctgag");
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);

    EXPECT_EQ(search_states.size(), 1);

    auto result = search_states.front().traversed_path;
    VariantSitePath expected = {
            VariantLocus {5, 2}
    };
    EXPECT_EQ(result, expected);
}


TEST(SearchStates, ReadCrossingTwoAlleles_CorrectVarLocusPaths) {
    Sequence kmer = encode_dna_bases("tct");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7T8c8CT", kmers);

    auto read = encode_dna_bases("cagtct");
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto traversed_path = search_states.front().traversed_path;
    VariantSitePath expected_traversed = { VariantLocus {7, 1}};
    EXPECT_EQ(traversed_path, expected_traversed);

    auto traversing_path = search_states.front().traversing_path;
    VariantSitePath expected_traversing = {VariantLocus{5, ALLELE_UNKNOWN}};
    EXPECT_EQ(traversing_path, expected_traversing);
}


TEST(SearchStates, StartWithinAlleleEndWithinAnother_CorrectVarLocusPath) {
    Sequence kmer = encode_dna_bases("gag");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7GAG8c8ct", kmers);

    auto read = encode_dna_bases("caggag");
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto traversed_path = search_states.front().traversed_path;
    VariantSitePath expected_traversed = { VariantLocus {7, 1}};
    EXPECT_EQ(traversed_path, expected_traversed);

    auto traversing_path = search_states.front().traversing_path;
    VariantSitePath expected_traversing = {VariantLocus{5, ALLELE_UNKNOWN}};
    EXPECT_EQ(traversing_path, expected_traversing);
}


/*
 * A case where we end the read mapping inside several alleles of the same site.
 * We test: correct indexing, correct base extension, correct allele id specification.
 */
TEST(MultiStep, RunIndexingExtensionIdSpecification_CorrectOutputs) {
    Sequence kmer = encode_dna_bases("c");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5gC6aC6C6t6Cg", kmers);

    // We expect five occurrences of 'C' at this stage, in a single SA interval
    auto search_states = setup.kmer_index.at(kmer);
    EXPECT_EQ(search_states.size(), 1);
    SA_Interval sa = search_states.front().sa_interval;
    EXPECT_EQ(sa.second - sa.first + 1, 5);

    // Next up, look for a C
    int_Base pattern_char = 2;
    search_states = process_read_char_search_states(pattern_char,
                                                    search_states,
                                                    setup.prg_info);

    // concurrent allele querying
    // Expect three occurrences of 'CC' at this stage, in a single SA interval
    EXPECT_EQ(search_states.size(), 1);
    EXPECT_EQ(search_states.front().traversing_path.back().second, ALLELE_UNKNOWN);
    EXPECT_EQ(search_states.front().sa_interval.second - search_states.front().sa_interval.first + 1, 3);
}

TEST(SearchStates, OneMappingEncapsulatedByAllele) {
    Sequence kmer = encode_dna_bases("tagt");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("t5c6gCTTAGT6aa", kmers);


    auto read = encode_dna_bases("cttagt");
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto result = search_states.front().variant_site_state;
    SearchVariantSiteState expected = SearchVariantSiteState::within_variant_site;
    EXPECT_EQ(result, expected);

    VariantLocus cov = {5, 2};
    EXPECT_EQ(search_states.front().traversed_path.front(), cov);
}

TEST(SearchStates, StartAndEndInSite_CorrectSearchStates) {
    Sequence kmer = encode_dna_bases("tagt");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup("t5c6gcttagtacgcttagt6aa", kmers);

    auto read = encode_dna_bases("cttagt");
    auto result = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);

    SearchStates expected = {
            SearchState{
                    SA_Interval{7, 8},
                    VariantSitePath{
                            VariantLocus{5, 2}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            }
    };

    EXPECT_EQ(result, expected);
}


TEST(SearchStates_Nested, MapIntoAndOutOfNestedSite_CorrectSearchStates) {
    Sequence kmer = encode_dna_bases("ac");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup_nested("a[c,g[ct,t]a]c", kmers);

    auto read = encode_dna_bases("agtac");
    auto result = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);

    SearchStates expected = {
            SearchState{
                SA_Interval{1, 1},
                VariantSitePath{
                    VariantLocus{7, 2},
                    VariantLocus{5, 2}
                },
                VariantSitePath{},
                SearchVariantSiteState::outside_variant_site
            }
    };
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
    Sequence kmer = encode_dna_bases("t");
    Sequences kmers = {kmer};
    prg_setup setup;
    setup.setup_nested("t[a[c,g][c,g],]t", kmers);

    auto read = encode_dna_bases("tt");
    auto result_direct_deletion = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);

    SearchStates expected_direct_deletion = {
            SearchState{
                SA_Interval{7, 7},
                VariantSitePath{VariantLocus{5, 2}},
                VariantSitePath{},
                SearchVariantSiteState::outside_variant_site
            }
    };
    EXPECT_EQ(result_direct_deletion, expected_direct_deletion);

    auto read2 = encode_dna_bases("tacct");
    auto result_exit_entry = search_read_backwards(read2, kmer, setup.kmer_index, setup.prg_info);

    SearchStates expected_exit_entry = {
           SearchState{
               SA_Interval{7, 7},
               VariantSitePath{
                   VariantLocus{9, 1},
                   VariantLocus{7, 1},
                   VariantLocus{5, 1}
               },
               VariantSitePath{},
               SearchVariantSiteState::outside_variant_site
           }
    };
    EXPECT_EQ(result_exit_entry, expected_exit_entry);
}

class Coverage_Nested_DoubleNesting : public ::testing::Test {
    // Double nesting meaning a bubble inside a bubble inside a bubble
protected:
    void SetUp() {
        Sequence kmer = encode_dna_bases("TA");
        Sequences kmers = {kmer};
        setup.setup_nested("A[[A[CCC,c],t],g]TA", kmers);
    }
    prg_setup setup;
    prg_positions positions{0, 3, 5, 9, 12, 15, 17}; // All the nodes in the cov graph with sequence

    Sequence read1 = encode_dna_bases("AACCCTA");
    Sequence read2 = encode_dna_bases("CTA");
};

TEST_F(Coverage_Nested_DoubleNesting, ReadEndsInsideNestedSite_CorrectCoverage){
    // PRG: "A[[A[CCC,c],t],g]TA"; Read: "AACCCTA"
    quasimap_read(read1, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
    // The read is compatible with the first allele of all three sites in the PRG
    SitesGroupedAlleleCounts expectedGpAlCounts = {
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
    };
    EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

    auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
    AlleleCoverage expectedPbCov{
      BaseCoverage{}, BaseCoverage{1}, BaseCoverage{1, 1, 1},
      BaseCoverage{0}, BaseCoverage{0}, BaseCoverage{0},
      BaseCoverage{}
    };
    EXPECT_EQ(PbCov, expectedPbCov);
}

TEST_F(Coverage_Nested_DoubleNesting, ReadMultiMaps_CorrectCoverage){
    // PRG: "A[[A[CCC,c],t],g]TA"; Read: "CTA"
    quasimap_read(read2, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
    // The read is compatible with the two alleles of the most nested site in the PRG string
    SitesGroupedAlleleCounts expectedGpAlCounts = {
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {{AlleleIds {0, 1}, 1}},
    };
    EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

    auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
    AlleleCoverage expectedPbCov{
            BaseCoverage{}, BaseCoverage{0}, BaseCoverage{0, 0, 1},
            BaseCoverage{1}, BaseCoverage{0}, BaseCoverage{0},
            BaseCoverage{}
    };
    EXPECT_EQ(PbCov, expectedPbCov);
}


class Coverage_Nested_SingleNestingPlusSNP : public ::testing::Test {
protected:
    void SetUp() {
        Sequences kmers = {
                encode_dna_bases("C"),
                encode_dna_bases("G"),
                encode_dna_bases("T"),
        };
        setup.setup_nested("a[t[tt,t]t,a[at,]a]g[c,g]", kmers);
    }
    prg_setup setup;
    prg_positions positions{0, 2, 4, 7, 9, 11, 13, 17, 19, 21, 23}; // All the nodes in the cov graph with sequence

    Sequence read1 = encode_dna_bases("ATTTTGC");
    Sequence read2 = encode_dna_bases("TT");
    Sequence read3 = encode_dna_bases("AAAGG");
};

TEST_F(Coverage_Nested_SingleNestingPlusSNP, FullyCrossingRead_CorrectCoverage){
    //PRG: "A[T[TT,T]T,a[at,]a]G[C,g]" ; Read: "ATTTTGC"
    quasimap_read(read1, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
    // The read is compatible with the two alleles of the most nested site in the PRG string
    SitesGroupedAlleleCounts expectedGpAlCounts = {
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {},
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
    };
    EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

    auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
    AlleleCoverage expectedPbCov{
            BaseCoverage{},
            BaseCoverage{1}, BaseCoverage{1, 1}, BaseCoverage{0}, BaseCoverage{1},
            BaseCoverage{0}, BaseCoverage{0, 0}, BaseCoverage{0},
            BaseCoverage{},
            BaseCoverage{1}, BaseCoverage{0}
    };
    EXPECT_EQ(PbCov, expectedPbCov);
}

TEST_F(Coverage_Nested_SingleNestingPlusSNP, VeryMultiMappingRead_CorrectCoverage) {
    //PRG: "A[T[TT,T]T,a[at,]a]G[C,g]" ; Read: "TT"
    // This read should have 5 mapping instances: one is encapsulated(=empty traversing and traversed),
    // two are in 'traversing' mode, two are in 'traversed' mode. All are encapsulated inside site 0 as well.

    quasimap_read(read2, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
    // The read is compatible with the two alleles of the most nested site in the PRG string
    SitesGroupedAlleleCounts expectedGpAlCounts = {
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {{AlleleIds {0, 1}, 1}},
            GroupedAlleleCounts {},
            GroupedAlleleCounts {},
    };
    EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

    auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
    AlleleCoverage expectedPbCov{
            BaseCoverage{},
            BaseCoverage{1}, BaseCoverage{1, 1}, BaseCoverage{1}, BaseCoverage{1},
            BaseCoverage{0}, BaseCoverage{0, 0}, BaseCoverage{0},
            BaseCoverage{},
            BaseCoverage{0}, BaseCoverage{0}
    };
    EXPECT_EQ(PbCov, expectedPbCov);
}

TEST_F(Coverage_Nested_SingleNestingPlusSNP, MapThroughDirectDeletion_CorrectCoverage) {
    //PRG: "A[t[tt,t]t,A[at,]A]G[c,G]" ; Read: "AAAGG"
    quasimap_read(read3, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &GpAlCounts = setup.coverage.grouped_allele_counts;
    SitesGroupedAlleleCounts expectedGpAlCounts = {
            GroupedAlleleCounts {{AlleleIds {1}, 1}},
            GroupedAlleleCounts {},
            GroupedAlleleCounts {{AlleleIds {1}, 1}},
            GroupedAlleleCounts {{AlleleIds {1}, 1}},
    };
    EXPECT_EQ(GpAlCounts, expectedGpAlCounts);

    auto PbCov = collect_coverage(setup.prg_info.coverage_graph, positions);
    AlleleCoverage expectedPbCov{
            BaseCoverage{},
            BaseCoverage{0}, BaseCoverage{0, 0}, BaseCoverage{0}, BaseCoverage{0},
            BaseCoverage{1}, BaseCoverage{0, 0}, BaseCoverage{1},
            BaseCoverage{},
            BaseCoverage{0}, BaseCoverage{1}
    };
    EXPECT_EQ(PbCov, expectedPbCov);
}
