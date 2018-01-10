#include <cctype>

#include "gtest/gtest.h"

#include "test_utils.hpp"
#include "kmer_index.hpp"
#include "coverage_analysis.hpp"


TEST(CoverageAnalysis, GivenOneVariantSite_CorrectAlleleCoverageStructure) {
    const auto prg_raw = "gcgct5gg6agtg5ctgt";
    const auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_allele_coverage_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, GivenTwoVariantSite_CorrectAlleleCoverageStructure) {
    const auto prg_raw = "gcgct5gg6agtg5cccc7t8g7t";
    const auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_allele_coverage_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, GivenThreeVariantSites_CorrectAlleleCoverageStructure) {
    const auto prg_raw = "5gg6agtg5c7t8g8c7t9ccccc10t9";
    const auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_allele_coverage_structure(prg_info);
    AlleleSumCoverage expected = {
            {0, 0},
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, GivenReadAndKmerSize_CorrectKmerReturned) {
    auto read = encode_dna_bases("accgaatt");
    uint32_t kmer_size = 3;
    auto result = get_kmer_from_read(kmer_size, read);
    auto expected = encode_dna_bases("att");
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("gccta");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("agccta");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("agtcta");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("ctgagtcta");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 0},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadCrossingMultipleVariantSitesEndingInAllele_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagtcta");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 1},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, NonMappingReadCrossingAllele_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tgtcta");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadEndsInAllele_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("ctc");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("gctc");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadStartsInAllele_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagt");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 1},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadWithNoMatchingKmer_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagc");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadMapsToThreePositions_CorrectAlleleCoverage) {
    const auto prg_raw = "tag5tc6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagt");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 2},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadEntierlyWithinAllele_CoverageNotRecorded) {
    const auto prg_raw = "gct5cccc6g6t5ag";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("ccc");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("cccc");

    quasimap_read(read,
                  coverage,
                  kmer_index,
                  prg_info,
                  parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, MappingMultipleIdenticalReads_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("tagt"),
            encode_dna_bases("tagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read,
                      coverage,
                      kmer_index,
                      prg_info,
                      parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 2},
            {2, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, MappingTwoReadsIdenticalKmers_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read,
                      coverage,
                      kmer_index,
                      prg_info,
                      parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 1, 1},
            {2, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, MappingThreeReadsIdenticalKmers_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Pattern kmer = encode_dna_bases("agt");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt"),
            encode_dna_bases("cagt")
    };

    for (const auto &read: reads) {
        quasimap_read(read,
                      coverage,
                      kmer_index,
                      prg_info,
                      parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 1, 1},
            {3, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, MappingThreeReadsDifferentKmers_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("agt"),
            encode_dna_bases("agc"),
    };
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("gagt"),
            encode_dna_bases("tagt"),
            encode_dna_bases("cagc")
    };

    for (const auto &read: reads) {
        quasimap_read(read,
                      coverage,
                      kmer_index,
                      prg_info,
                      parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 1, 1},
            {2, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, MappingThreeReadsOneReadMappsTwice_CorrectAlleleCoverage) {
    const auto prg_raw = "gcac5t6g6c5ta7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto coverage = generate_coverage_info_structure(prg_info);

    Patterns kmers = {
            encode_dna_bases("cta"),
            encode_dna_bases("act"),
    };
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("accta"),
            encode_dna_bases("gcact"),
    };

    for (const auto &read: reads) {
        quasimap_read(read,
                      coverage,
                      kmer_index,
                      prg_info,
                      parameters);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1},
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}