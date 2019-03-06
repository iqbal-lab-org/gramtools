#include <cctype>

#include "gtest/gtest.h"

#include "../test_utils.hpp"
#include "kmer_index/build.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"


using namespace gram;


TEST(Quasimap, GivenReadAndKmerSize_CorrectKmerReturned) {
    auto read = encode_dna_bases("accgaatt");
    uint32_t kmer_size = 3;
    auto result = get_kmer_from_read(kmer_size, read);
    auto expected = encode_dna_bases("att");
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "tag5tc6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5cccc6g6t5ag";
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
PRG: ac5t6cagtagtc5ta
i	F	BWT	text	SA	suffix
0	0	1	1	    16	0
1	1	4	2	    15	1 0
2	1	0	5	    0	1 2 5 4 6 2 1 3 4 1 3 4 2 5 4 1 0
3	1	2	4	    6	1 3 4 1 3 4 2 5 4 1 0
4	1	4	6	    9	1 3 4 2 5 4 1 0
5	2	6	2	    5	2 1 3 4 1 3 4 2 5 4 1 0
6	2	4	1	    12	2 5 4 1 0
7	2	1	3	    1	2 5 4 6 2 1 3 4 1 3 4 2 5 4 1 0
8	3	1	4	    7	3 4 1 3 4 2 5 4 1 0
9	3	1	1	    10	3 4 2 5 4 1 0
10	4	5	3	    14	4 1 0
11	4	3	4	    8	4 1 3 4 2 5 4 1 0
12	4	3	2	    11	4 2 5 4 1 0
13	4	5	5	    3	4 6 2 1 3 4 1 3 4 2 5 4 1 0
14	5	2	4	    13	5 4 1 0
15	5	2	1	    2	5 4 6 2 1 3 4 1 3 4 2 5 4 1 0
16	6	4	0	    4	6 2 1 3 4 1 3 4 2 5 4 1 0
*/
TEST(Quasimap, ReadMapsWithinAllele_SumCoverageIsOne) {
    auto prg_raw = "ac5t6cagtagtc5ta";
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
    auto prg_raw = "ac5t6cagtagttttgtagtc5ta";
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
    auto prg_raw = "gtagtac5gtagtact6t5ta";
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
    auto prg_raw = "tac5gta6gtt5ta";
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
    auto prg_raw = "c5ccc6agt6ccgt5taa";
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
    auto prg_raw = "ac5gtagtact6t6gggtagt5ta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
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
    auto prg_raw = "gcac5t6g6c5ta7t8c7cta";
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
