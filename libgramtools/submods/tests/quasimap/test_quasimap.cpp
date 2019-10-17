/**
 * @file
 * Test high-level quasimapping routine: searching for full kmers or full reads
 * The tested outputs are recorded `Coverage`s and `SearchState` internals (SA intervals, `VariantLocus` paths)
 */

#include "gtest/gtest.h"

#include "src_common/generate_prg.hpp"
#include "kmer_index/build.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"
#include "common/utils.hpp"


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
        size_t kmer_size = kmers.front().size();
       for (auto const& kmer : kmers) assert(kmer_size == kmer.size());

       auto encoded_prg = encode_prg(raw_prg);
        prg_info = generate_prg_info(encoded_prg);
        // TODO: the calls to rank_support setup in `generate_prg_info` do not set up the rank support properly
        //  when leaving this scope.
        prg_info.rank_bwt_a = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_a);
        prg_info.rank_bwt_c = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_c);
        prg_info.rank_bwt_g = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_g);
        prg_info.rank_bwt_t = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_t);
        prg_info.prg_markers_rank = sdsl::rank_support_v<1>(&prg_info.prg_markers_mask);
        prg_info.prg_markers_select = sdsl::select_support_mcl<1>(&prg_info.prg_markers_mask);

        coverage = coverage::generate::empty_structure(prg_info);

        parameters.kmers_size = kmer_size;
        kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);
    }
};

TEST(ReadQuasimap, GivenReadAndKmerSize_CorrectKmerReturned) {
    auto read = encode_dna_bases("accgaatt");
    uint32_t kmer_size = 3;
    auto result = get_kmer_from_read(kmer_size, read);
    auto expected = encode_dna_bases("att");
    EXPECT_EQ(result, expected);
}


TEST(ReadQuasimap, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, ReadCrossTwoSitesAndEndsInSite_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, ReadDoesNotMap_EmptyAlleleCoverage) {
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


TEST(ReadQuasimap, ReadEndsInAllele_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, ReadStartsInAllele_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, ReadWithNoMatchingKmer_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, ReadMapsToThreePositions_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, ReadEntierlyWithinAllele_CoverageRecorded) {
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

TEST(ReadQuasimap, ReadMapsWithinAllele_SumCoverageIsOne) {
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


TEST(ReadQuasimap, ReadMapsTwiceWithinAllele_SumCoverageIsOne) {
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


TEST(ReadQuasimap, ReadMapsWithinAlleleAndOutsideSite_CorrectSumCoverage) {
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


TEST(ReadQuasimap, ReadEndWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
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


TEST(ReadQuasimap, ReadStartWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
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


TEST(ReadQuasimap, EncapsulatedWithinTwoDifferentAlleles_CorrectAlleleSumCoverage) {
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


TEST(ReadQuasimap, MappingMultipleIdenticalReads_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, MappingTwoReadsIdenticalKmers_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, MappingThreeReadsIdenticalKmers_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, MappingThreeReadsDifferentKmers_CorrectAlleleCoverage) {
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


TEST(ReadQuasimap, MappingThreeReadsOneReadMappsTwice_CorrectAlleleCoverage) {
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

TEST(ReadQuasimap, StartOufofSiteAndEndInSite_CorrectSearchState) {
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

TEST(ReadQuasimap, StartInSiteAndMapOut_CorrectVarLocusPath) {
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


TEST(ReadQuasimap, StartOutOfSiteAndMapThrough_CorrectVarLocusPath) {
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


TEST(ReadQuasimap, ReadCrossingTwoAlleles_CorrectVarLocusPath) {
    Pattern kmer = encode_dna_bases("tct");
    Patterns kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7T8c8CT", kmers);

    auto read = encode_dna_bases("cagtct");
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto result = search_states.front().traversed_path;
    VariantSitePath expected = {
            VariantLocus {7, 1},
            VariantLocus {5, 1},
    };
    EXPECT_EQ(result, expected);
}


TEST(ReadQuasimap, StartWithinAlleleEndWithinAnother_CorrectVarLocusPath) {
    Pattern kmer = encode_dna_bases("gag");
    Patterns kmers = {kmer};
    prg_setup setup;
    setup.setup("gct5c6g6t6ag7GAG8c8ct", kmers);

    auto read = encode_dna_bases("caggag");
    auto search_states = search_read_backwards(read, kmer, setup.kmer_index, setup.prg_info);
    EXPECT_EQ(search_states.size(), 1);

    auto result = search_states.front().traversed_path;
    VariantSitePath expected = {
            VariantLocus {7, 1},
            VariantLocus {5, 1},
    };
    EXPECT_EQ(result, expected);
}


/*
 * A case where we end the read mapping inside several alleles of the same site.
 * We test: correct indexing, correct base extension, correct allele id specification.
 */
TEST(MultiStepQuasimap, RunIndexingExtensionIdSpecification_CorrectOutputs) {
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
    // Expect three occurrences of 'CC' at this stage, in a single SA interval - because
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

TEST(ReadQuasimap, OneMappingEncapsulatedByAllele) {
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

TEST(ReadQuasimap, StartAndEndInSite_CorrectSearchStates) {
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

