/**
 * @file
 * Test high-level quasimapping routine: searching for full kmers or full reads.
 * Assessing results is in terms of SearchStates produced or coverage recorded.
 *
 * Suites:
 *  - SearchStates: test that you produce the right search states
 *  - AlleleSum: test that mapping increments the right allele sum coverage
 *  - GpedAlCounts: test that mapping increments the right grouped allele counts coverage
 *
 *  A "_Nested" prefix is added for nested PRGs.
 *
 */

#include "gtest/gtest.h"

#include "src_common/generate_prg.hpp"
#include "kmer_index/build.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"
#include "common/utils.hpp"
#include "quasimap/search/BWT_search.hpp"

using namespace gram;

class prg_setup{
public:
    PRG_Info prg_info;
    Coverage coverage;
    Parameters parameters;
    KmerIndex kmer_index;

    explicit prg_setup() {};
    void setup(std::string raw_prg,
            Patterns kmers){
       auto encoded_prg = encode_prg(raw_prg);
       internal_setup(encoded_prg, kmers);
    }

    void setup_nested(std::string raw_prg,
               Patterns kmers){
        auto encoded_prg = prg_string_to_ints(raw_prg);
        internal_setup(encoded_prg, kmers);
    }

private:
    void internal_setup(marker_vec encoded_prg, Patterns kmers){
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

TEST(GetKmer, GivenReadAndKmerSize_CorrectKmerReturned) {
    auto read = encode_dna_bases("accgaatt");
    uint32_t kmer_size = 3;
    auto result = get_kmer_from_read(kmer_size, read);
    auto expected = encode_dna_bases("att");
    EXPECT_EQ(result, expected);
}


TEST(AlleleSum, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("gccta");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadCrossTwoSitesAndEndsInSite_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadDoesNotMap_EmptyAlleleCoverage) {
    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadEndsInAllele_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("ctc");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadStartsInAllele_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadWithNoMatchingKmer_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadMapsToThreePositions_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
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


TEST(AlleleSum, ReadEntierlyWithinAllele_CoverageRecorded) {
    Pattern kmer = encode_dna_bases("ccc");
    Patterns kmers = {kmer};
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


/*
PRG: AC5T6CAGTAGTC6TA
i	BWT	SA	text_suffix
0	A	16
1	T	15	A
2	0	0	A C 5 T 6 C A G T A G T C 6 T A
3	C	6	A G T A G T C 6 T A
4	T	9	A G T C 6 T A
5	6	5	C A G T A G T C 6 T A
6	A	1	C 5 T 6 C A G T A G T C 6 T A
7	T	12	C 6 T A
8	A	7	G T A G T C 6 T A
9	A	10	G T C 6 T A
10	6	14	T A
11	G	8	T A G T C 6 T A
12	G	11	T C 6 T A
13	5	3	T 6 C A G T A G T C 6 T A
14	C	2	5 T 6 C A G T A G T C 6 T A
15	T	4	6 C A G T A G T C 6 T A
16	C	13	6 T A
*/

TEST(AlleleSum, ReadMapsWithinAllele_SumCoverageIsOne) {
    Patterns kmers = {
            encode_dna_bases("agt"),
    };
    prg_setup setup;
    setup.setup("ac5t6cagtagtc6ta",kmers);

    Pattern read = encode_dna_bases("gtagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSum, ReadMapsTwiceWithinAllele_SumCoverageIsOne) {
    Patterns kmers = {
            encode_dna_bases("agt"),
    };
    prg_setup setup;
    setup.setup("ac5t6cagtagttttgtagtc6ta",kmers);
    setup.parameters.seed = 42;

    Pattern read = encode_dna_bases("gtagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSum, ReadMapsWithinAlleleAndOutsideSite_CorrectSumCoverage) {
    Patterns kmers = {
            encode_dna_bases("agt"),
    };
    prg_setup setup;
    setup.setup("gtagtac5gtagtact6t6ta",kmers);
    setup.parameters.seed = 39;

    Pattern read = encode_dna_bases("gtagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSum, ReadEndWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
    Patterns kmers = {
            encode_dna_bases("cgt"),
    };
    prg_setup setup;
    setup.setup("tac5gta6gtt6ta", kmers);

    Pattern read = encode_dna_bases("tacgt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSum, ReadStartWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
    Patterns kmers = {
            encode_dna_bases("taa"),
    };
    prg_setup setup;
    setup.setup("c5ccc6agt6ccgt6taa", kmers);
    setup.parameters.seed = 39;

    Pattern read = encode_dna_bases("gttaa");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSum, EncapsulatedWithinTwoDifferentAlleles_CorrectAlleleSumCoverage) {
    Patterns kmers = {
            encode_dna_bases("agt"),
    };
    prg_setup setup;
    setup.setup("ac5gtagtact6t6gggtagt6ta", kmers);
    setup.parameters.seed = 42;

    Pattern read = encode_dna_bases("gtagt");
    quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSum, MappingMultipleIdenticalReads_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta", kmers);

    Patterns reads = {
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
}


TEST(AlleleSum, MappingTwoReadsIdenticalKmers_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta", kmers);

    Patterns reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);
    }

    const auto &result = setup.coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 1},
            {2, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(AlleleSum, MappingThreeReadsIdenticalKmers_CorrectAlleleCoverage) {
    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta", kmers);

    Patterns reads = {
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
}


TEST(AlleleSum, MappingThreeReadsDifferentKmers_CorrectAlleleCoverage) {
    Patterns kmers = {
            encode_dna_bases("agt"),
            encode_dna_bases("agc"),
    };
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7t8c8cta", kmers);

    Patterns reads = {
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


TEST(AlleleSum, MappingThreeReadsOneReadMappsTwice_CorrectAlleleCoverage) {
    Patterns kmers = {
            encode_dna_bases("cta"),
            encode_dna_bases("act"),
    };
    prg_setup setup;
    setup.setup("gcac5t6g6c6ta7t8c8cta", kmers);
    setup.parameters.seed = 42;

    Patterns reads = {
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

    auto read = encode_dna_bases("tagtaa");
    Pattern kmer = encode_dna_bases("gtaa");
    Patterns kmers = {kmer};
    auto kmer_size = 4;
    auto kmer_index = index_kmers(kmers, kmer_size, prg_info);

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
    auto search_state = final_search_states.front();
    const auto &result = search_state.traversed_path;
    VariantSitePath expected = {
            VariantLocus{5, 2},
    };
    EXPECT_EQ(result, expected);
}

TEST(SearchStates, StartOufofSiteAndEndInSite_CorrectSearchState) {
    Pattern kmer = encode_dna_bases("gtcc");
    Patterns kmers = {kmer};
    prg_setup setup;
    setup.setup("gcgct5c6g6T6AGTCCt", kmers);

    auto read = encode_dna_bases("tagtcc");
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
    EXPECT_EQ(search_states.size(), 1);

    //Do we end up in right place in SA index?
    auto search_state = search_states.front();
    auto result = search_state.sa_interval;
    SA_Interval expected = {14, 14};
    EXPECT_EQ(result, expected);

    auto path_result = search_state.traversed_path;
    VariantSitePath path_expected = {
            VariantLocus {5, 3} // We expect it to be traversed because we fully mapped the read, so sites got assigned.
    };
    EXPECT_EQ(path_result, path_expected);
}

TEST(SearchStates, StartInSiteAndMapOut_CorrectVarLocusPath) {
    Pattern kmer = encode_dna_bases("gctc");
    Patterns kmers = {kmer};
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
    Pattern kmer = encode_dna_bases("tgag");
    Patterns kmers = {kmer};
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
    Pattern kmer = encode_dna_bases("tct");
    Patterns kmers = {kmer};
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
    Pattern kmer = encode_dna_bases("gag");
    Patterns kmers = {kmer};
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
    Pattern kmer = encode_dna_bases("c");
    Patterns kmers = {kmer};
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

    // allele id specification
    // we should now have three search states of SA Interval size 1, each with a different traversed allele id
    gram::set_allele_ids(search_states, setup.prg_info);
    EXPECT_EQ(search_states.size(), 3);

    std::set<AlleleId> ids;
    for (auto const& search_state: search_states){
        SA_Interval sa = search_state.sa_interval;
        EXPECT_EQ(sa.second - sa.first + 1, 1);
        ids.insert(search_state.traversed_path.back().second);
    }
    std::set<AlleleId> expected{1,2,3};
    EXPECT_EQ(ids, expected);
}

TEST(SearchStates, OneMappingEncapsulatedByAllele) {
    Pattern kmer = encode_dna_bases("tagt");
    Patterns kmers = {kmer};
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
    Pattern kmer = encode_dna_bases("tagt");
    Patterns kmers = {kmer};
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
    Pattern kmer = encode_dna_bases("ac");
    Patterns kmers = {kmer};
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
    Pattern kmer = encode_dna_bases("t");
    Patterns kmers = {kmer};
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

TEST(GpedAlCounts_Nested, DoubleNesting_CorrectCoverage){
    Pattern kmer = encode_dna_bases("CTA");
    Patterns kmers = {kmer};
    prg_setup setup;
    setup.setup_nested("A[[A[CCC,c],t],g]TA", kmers);

    auto read1 = encode_dna_bases("AACCCTA");
    quasimap_read(read1, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);

    const auto &result = setup.coverage.grouped_allele_counts;
    // The read is compatible with the first allele of all three sites in the PRG
    SitesGroupedAlleleCounts expected = {
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
    };
    EXPECT_EQ(result, expected);

    // This read is also compatible with the same sites above
    auto read2 = encode_dna_bases("CCTA");

    setup.coverage = coverage::generate::empty_structure(setup.prg_info); // Clear out coverage
    EXPECT_NE(result, expected); // Make sure coverage has been invalidated

    quasimap_read(read2, setup.coverage, setup.kmer_index, setup.prg_info, setup.parameters);
    EXPECT_EQ(result, expected);
}
