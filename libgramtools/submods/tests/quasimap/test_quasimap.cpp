#include <cctype>

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

TEST(Quasimap, GivenReadAndKmerSize_CorrectKmerReturned) {
    auto read = encode_dna_bases("accgaatt");
    uint32_t kmer_size = 3;
    auto result = get_kmer_from_read(kmer_size, read);
    auto expected = encode_dna_bases("att");
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
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


TEST(Quasimap, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
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


TEST(Quasimap, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
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


TEST(Quasimap, ReadCrossTwoSitesAndEndsInSite_CorrectAlleleCoverage) {
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


TEST(Quasimap, ReadDoesNotMap_EmptyAlleleCoverage) {
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


TEST(Quasimap, ReadEndsInAllele_CorrectAlleleCoverage) {
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


TEST(Quasimap, ReadStartsInAllele_CorrectAlleleCoverage) {
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


TEST(Quasimap, ReadWithNoMatchingKmer_CorrectAlleleCoverage) {
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


TEST(Quasimap, ReadMapsToThreePositions_CorrectAlleleCoverage) {
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


TEST(Quasimap, ReadEntierlyWithinAllele_CoverageRecorded) {
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

TEST(Quasimap, ReadMapsWithinAllele_SumCoverageIsOne) {
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


TEST(Quasimap, ReadMapsTwiceWithinAllele_SumCoverageIsOne) {
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


TEST(Quasimap, ReadMapsWithinAlleleAndOutsideSite_CorrectSumCoverage) {
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


TEST(Quasimap, ReadEndWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
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


TEST(Quasimap, ReadStartWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
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


TEST(Quasimap, EncapsulatedWithinTwoDifferentAlleles_CorrectAlleleSumCoverage) {
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


TEST(Quasimap, MappingMultipleIdenticalReads_CorrectAlleleCoverage) {
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


TEST(Quasimap, MappingTwoReadsIdenticalKmers_CorrectAlleleCoverage) {
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


TEST(Quasimap, MappingThreeReadsIdenticalKmers_CorrectAlleleCoverage) {
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


TEST(Quasimap, MappingThreeReadsDifferentKmers_CorrectAlleleCoverage) {
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


TEST(Quasimap, MappingThreeReadsOneReadMappsTwice_CorrectAlleleCoverage) {
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

