#include <cctype>

#include "gtest/gtest.h"

#include "src_common/generate_prg.hpp"
#include "kmer_index/build.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"
#include "common/utils.hpp"


using namespace gram;


TEST(Quasimap, GivenReadAndKmerSize_CorrectKmerReturned) {
    auto read = encode_dna_bases("accgaatt");
    uint32_t kmer_size = 3;
    auto result = get_kmer_from_read(kmer_size, read);
    auto expected = encode_dna_bases("att");
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6aG7t8C8CTA";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("gccta");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("agccta");

    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("agtcta");

    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("ctgagtcta");

    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 0},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadCrossingMultipleVariantSitesEndingInAllele_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagtcta");

    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 1},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, NonMappingReadCrossingAllele_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tgtcta");

    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadEndsInAllele_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("ctc");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("gctc");

    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadStartsInAllele_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagt");

    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 1},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadWithNoMatchingKmer_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagc");

    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadMapsToThreePositions_CorrectAlleleCoverage) {
    auto prg_raw = "tag5tc6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 3;
    parameters.seed = 42;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagt");
    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadEntierlyWithinAllele_CoverageRecorded) {
    auto prg_raw = "gct5cccc6g6t6ag";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("ccc");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("cccc");
    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
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
    auto prg_raw = "ac5t6cagtagtc6ta";
    auto prg_info = generate_prg_info(prg_raw);

    auto coverage = coverage::generate::empty_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("agt"),
    };
    Parameters parameters = {};
    parameters.kmers_size = 3;
    parameters.seed = 42;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Pattern read = encode_dna_bases("gtagt");
    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadMapsTwiceWithinAllele_SumCoverageIsOne) {
    auto prg_raw = "ac5t6cagtagttttgtagtc6ta";
    auto prg_info = generate_prg_info(prg_raw);

    auto coverage = coverage::generate::empty_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("agt"),
    };
    Parameters parameters = {};
    parameters.kmers_size = 3;
    parameters.seed = 42;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Pattern read = encode_dna_bases("gtagt");
    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadMapsWithinAlleleAndOutsideSite_CorrectSumCoverage) {
    auto prg_raw = "gtagtac5gtagtact6t6ta";
    auto prg_info = generate_prg_info(prg_raw);

    auto coverage = coverage::generate::empty_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("agt"),
    };
    Parameters parameters = {};
    parameters.kmers_size = 3;
    parameters.seed = 39;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Pattern read = encode_dna_bases("gtagt");
    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadEndWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
    auto prg_raw = "tac5gta6gtt6ta";
    auto prg_info = generate_prg_info(prg_raw);

    auto coverage = coverage::generate::empty_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("cgt"),
    };
    Parameters parameters = {};
    parameters.kmers_size = 3;
    parameters.seed = 39;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Pattern read = encode_dna_bases("tacgt");
    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadStartWithinSingleSiteTwoAlleles_BothAlleleCoverage) {
    auto prg_raw = "c5ccc6agt6ccgt6taa";
    auto prg_info = generate_prg_info(prg_raw);

    auto coverage = coverage::generate::empty_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("taa"),
    };
    Parameters parameters = {};
    parameters.kmers_size = 3;
    parameters.seed = 39;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Pattern read = encode_dna_bases("gttaa");
    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, EncapsulatedWithinTwoDifferentAlleles_CorrectAlleleSumCoverage) {
    auto prg_raw = "ac5gtagtact6t6gggtagt6ta";
    auto prg_info = generate_prg_info(prg_raw);

    auto coverage = coverage::generate::empty_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("agt"),
    };
    Parameters parameters = {};
    parameters.kmers_size = 3;
    parameters.seed = 42;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Pattern read = encode_dna_bases("gtagt");
    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, MappingMultipleIdenticalReads_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("tagt"),
            encode_dna_bases("tagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read, coverage, kmer_index, prg_info, parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 2},
            {2, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, MappingTwoReadsIdenticalKmers_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read, coverage, kmer_index, prg_info, parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 1},
            {2, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, MappingThreeReadsIdenticalKmers_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters = {};
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt"),
            encode_dna_bases("cagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read, coverage, kmer_index, prg_info, parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 1, 1},
            {3, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, MappingThreeReadsDifferentKmers_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t6ag7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("agt"),
            encode_dna_bases("agc"),
    };
    Parameters parameters = {};
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt"),
            encode_dna_bases("cagc")
    };

    for (const auto &read: reads) {
        quasimap_read(read, coverage, kmer_index, prg_info, parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 1, 1},
            {2, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, MappingThreeReadsOneReadMappsTwice_CorrectAlleleCoverage) {
    auto prg_raw = "gcac5t6g6c6ta7t8c8cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("cta"),
            encode_dna_bases("act"),
    };
    Parameters parameters = {};
    parameters.kmers_size = 3;
    parameters.seed = 42;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("accta"),
            encode_dna_bases("gcact"),
    };

    for (const auto &read: reads) {
        quasimap_read(read, coverage, kmer_index,
                      prg_info, parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}

