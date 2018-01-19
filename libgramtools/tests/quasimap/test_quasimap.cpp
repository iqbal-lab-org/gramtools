#include <cctype>

#include "gtest/gtest.h"

#include "../test_utils.hpp"
#include "kmer_index/kmer_index.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"


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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagt");
    uint32_t random_seed = 42;
    quasimap_read(read, coverage, kmer_index, prg_info, parameters, random_seed);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, ReadEntierlyWithinAllele_CoverageNotRecorded) {
    auto prg_raw = "gct5cccc6g6t5ag";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    Pattern kmer = encode_dna_bases("ccc");
    Patterns kmers = {kmer};
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    const auto read = encode_dna_bases("cccc");
    quasimap_read(read, coverage, kmer_index, prg_info, parameters);

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {0, 0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(Quasimap, MappingMultipleIdenticalReads_CorrectAlleleCoverage) {
    auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
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
    Parameters parameters;
    parameters.kmers_size = 3;
    auto kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);

    Patterns reads = {
            encode_dna_bases("accta"),
            encode_dna_bases("gcact"),
    };

    uint32_t random_seed = 42;
    for (const auto &read: reads) {
        quasimap_read(read, coverage, kmer_index,
                      prg_info, parameters, random_seed);
    }

    const auto &result = coverage.allele_sum_coverage;
    AlleleSumCoverage expected = {
            {1, 0, 1},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}