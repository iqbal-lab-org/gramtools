#include <cctype>

#include "gtest/gtest.h"

#include "test_utils.hpp"
#include "kmer_index.hpp"
#include "coverage_analysis.hpp"


TEST(CoverageAnalysis, GivenOneVariantSite_CorrectAlleleCoverageStructure) {
    const auto prg_raw = "gcgct5gg6agtg5ctgt";
    const auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_allele_coverage_structure(prg_info);
    AlleleCoverage expected = {
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, GivenTwoVariantSite_CorrectAlleleCoverageStructure) {
    const auto prg_raw = "gcgct5gg6agtg5cccc7t8g7t";
    const auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_allele_coverage_structure(prg_info);
    AlleleCoverage expected = {
            {0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, GivenThreeVariantSites_CorrectAlleleCoverageStructure) {
    const auto prg_raw = "5gg6agtg5c7t8g8c7t9ccccc10t9";
    const auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_allele_coverage_structure(prg_info);
    AlleleCoverage expected = {
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
    auto allele_coverage = generate_allele_coverage_structure(prg_info);

    Pattern kmer = encode_dna_bases("gccta");
    Patterns kmers = {kmer};
    Parameters params;
    params.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, params.kmers_size, prg_info);

    const auto read = encode_dna_bases("agccta");

    quasimap_read(read,
                  allele_coverage,
                  kmer_index,
                  prg_info,
                  params);

    const auto &result = allele_coverage;
    AlleleCoverage expected = {
            {0, 0, 0},
            {0, 1}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto allele_coverage = generate_allele_coverage_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters params;
    params.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, params.kmers_size, prg_info);

    const auto read = encode_dna_bases("agtcta");

    quasimap_read(read,
                  allele_coverage,
                  kmer_index,
                  prg_info,
                  params);

    const auto &result = allele_coverage;
    AlleleCoverage expected = {
            {0, 0, 0},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST(CoverageAnalysis, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
    const auto prg_raw = "gct5c6g6t5ag7t8c7cta";
    const auto prg_info = generate_prg_info(prg_raw);
    auto allele_coverage = generate_allele_coverage_structure(prg_info);

    Pattern kmer = encode_dna_bases("gtcta");
    Patterns kmers = {kmer};
    Parameters params;
    params.kmers_size = 5;
    auto kmer_index = index_kmers(kmers, params.kmers_size, prg_info);

    const auto read = encode_dna_bases("tagtcta");

    quasimap_read(read,
                  allele_coverage,
                  kmer_index,
                  prg_info,
                  params);

    const auto &result = allele_coverage;
    AlleleCoverage expected = {
            {0, 0, 1},
            {1, 0}
    };
    EXPECT_EQ(result, expected);
}
