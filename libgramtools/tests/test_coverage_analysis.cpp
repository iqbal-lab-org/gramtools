#include <cctype>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "gtest/gtest.h"

#include "kmers.hpp"
#include "coverage_analysis.hpp"


class CoverageAnalysis : public ::testing::Test {

protected:
    std::string prg_fpath;
    std::string allele_coverage_fpath;

    void SetUp() override {
        boost::uuids::uuid uuid = boost::uuids::random_generator()();
        const auto uuid_str = boost::lexical_cast<std::string>(uuid);
        prg_fpath = "./prg_" + uuid_str;
        allele_coverage_fpath = "./allele_coverage_" + uuid_str;
    }

    void TearDown() override {
        std::remove(prg_fpath.c_str());
        std::remove(allele_coverage_fpath.c_str());
    }

    FM_Index fm_index_from_raw_prg(const std::string &prg_raw) {
        std::vector<uint64_t> prg = encode_prg(prg_raw);
        dump_encoded_prg(prg, prg_fpath);
        FM_Index fm_index;
        // TODO: constructing from memory with sdsl::construct_im appends 0 which corrupts
        sdsl::construct(fm_index, prg_fpath, 8);
        return fm_index;
    }

    PRG_Info generate_prg_info(const std::string &prg_raw) {
        PRG_Info prg_info;
        prg_info.fm_index = fm_index_from_raw_prg(prg_raw);
        prg_info.sites_mask = generate_sites_mask(prg_raw);
        // prg_info.dna_rank = calculate_ranks(prg_info.fm_index);
        prg_info.allele_mask = generate_allele_mask(prg_raw);
        prg_info.max_alphabet_num = max_alphabet_num(prg_raw);
        return prg_info;
    }

};


TEST_F(CoverageAnalysis, GivenOneVariantSite_CorrectAlleleCoverageStructure) {
    const auto prg_raw = "gcgct5gg6agtg5ctgt";
    const auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_allele_coverage_structure(prg_info);
    AlleleCoverage expected = {
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST_F(CoverageAnalysis, GivenTwoVariantSite_CorrectAlleleCoverageStructure) {
    const auto prg_raw = "gcgct5gg6agtg5cccc7t8g7t";
    const auto prg_info = generate_prg_info(prg_raw);

    auto result = generate_allele_coverage_structure(prg_info);
    AlleleCoverage expected = {
            {0, 0},
            {0, 0}
    };
    EXPECT_EQ(result, expected);
}


TEST_F(CoverageAnalysis, GivenThreeVariantSites_CorrectAlleleCoverageStructure) {
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


TEST_F(CoverageAnalysis, GivenReadAndKmerSize_CorrectKmerReturned) {
    auto read = encode_dna_bases("accgaatt");
    uint32_t kmer_size = 3;
    auto result = get_kmer_from_read(kmer_size, read);
    auto expected = encode_dna_bases("att");
    EXPECT_EQ(result, expected);
}


TEST_F(CoverageAnalysis, ReadCrossingSecondVariantSecondAllele_CorrectAlleleCoverage) {
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


TEST_F(CoverageAnalysis, ReadCrossingSecondVariantFirstAllele_CorrectAlleleCoverage) {
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


TEST_F(CoverageAnalysis, ReadCrossingMultipleVariantSites_CorrectAlleleCoverage) {
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


TEST_F(CoverageAnalysis, GivenAlleleCoverage_CorrectlyWrittenToFile) {
    AlleleCoverage allele_coverage = {
            {0, 0, 1},
            {1, 0}
    };
    Parameters params;
    params.allele_coverage_fpath = allele_coverage_fpath;
    dump_allele_coverage(allele_coverage, params);

    std::vector<std::string> result;
    std::string line;
    std::ifstream in_file_handle(allele_coverage_fpath);
    while (std::getline(in_file_handle, line)) {
        result.push_back(line);
    }
    std::vector<std::string> expected = {
            "0 0 1",
            "1 0"
    };
    EXPECT_EQ(result, expected);
}
