#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <sdsl/suffix_arrays.hpp>

#include "gtest/gtest.h"

#include "map.hpp"
#include "bidir_search_bwd.hpp"


std::vector<uint8_t> encode_read(const std::string &read) {
    std::vector<uint8_t> encoded_read;
    for (uint16_t i = 0; i < read.length(); i++) {
        if (read[i] == 'A' or read[i] == 'a') encoded_read.push_back(1);
        if (read[i] == 'C' or read[i] == 'c') encoded_read.push_back(2);
        if (read[i] == 'G' or read[i] == 'g') encoded_read.push_back(3);
        if (read[i] == 'T' or read[i] == 't') encoded_read.push_back(4);
    }
    return encoded_read;
}


class BidirSearchBackward : public ::testing::Test {

protected:
    std::string prg_fpath;

    virtual void SetUp() {
        boost::uuids::uuid uuid = boost::uuids::random_generator()();
        const auto uuid_str = boost::lexical_cast<std::string>(uuid);
        prg_fpath = "./prg_" + uuid_str;
    }

    virtual void TearDown() {
        std::remove(prg_fpath.c_str());
    }

    FM_Index fm_index_from_raw_prg(const std::string &prg_raw) {
        std::vector<uint64_t> prg = encode_prg(prg_raw);
        dump_encoded_prg(prg, prg_fpath);
        FM_Index fm_index;
        sdsl::construct(fm_index, prg_fpath, 8);
        return fm_index;
    }
};


TEST_F(BidirSearchBackward, MatchSingleCharecter) {
    const std::string prg_raw = "a";
    const uint64_t max_alphabet_num = 4;
    const std::string read = "a";
    const std::vector<int> allele_mask = {0};

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank &rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated, rank_all, fm_index);

    EXPECT_FALSE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {1, 2}
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected_sites = {
            {}
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchSingleVariantSiteOnly) {
    // aligns across SNP allele 1 (and both flanks)
    const std::string prg_raw = "catttacaca5g6t5aactagagagca";
    const uint64_t max_alphabet_num = 6;
    const std::string read = "ttacacagaactagagag";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 2, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank &rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated, rank_all, fm_index);

    EXPECT_TRUE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {22, 23}
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected = {
            {VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BidirSearchBackward, MatchTwoVariantSitesOnly) {
    const std::string prg_raw = "catttacaca5g6t5aactag7a8g7agcagggt";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "ttacacagaactagaagcag";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 2, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 2, 0, 0, 0, 0, 0,
            0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_TRUE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {26, 27}
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected = {
            {VariantSite(7, {1}), VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BidirSearchBackward, MatchTwoVariantSitesOnly_TwoVariantSitesIdentified) {
    const std::string prg_raw = "catttacaca5g6t5aactag7a8g7agcagggt";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "ttacacagaactagaagcag";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 2, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 2, 0, 0, 0, 0, 0,
            0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    const Sites expected = {
            {VariantSite(7, {1}), VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BidirSearchBackward, MatchTwoVariantSitesOnly_DeleteFirstIntervalTrue) {
    const std::string prg_raw = "catttacaca5g6t5aactag7a8g7agcagggt";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "ttacacagaactagaagcag";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 2, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 2, 0, 0, 0, 0, 0,
            0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_TRUE(delete_first_interval);
}


TEST_F(BidirSearchBackward, MatchOneVariantSiteMatchOneNonVariantSite) {
    //one match crosses allele 1, and the other in nonvar
    const std::string prg_raw = "catttacaca5g6t5aactagagagcaacagaactctct";
    const uint64_t max_alphabet_num = 6;
    const std::string read = "acagaac";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_FALSE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {5, 6},
            {6, 7},
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected_sites = {
            {},
            {VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchOneNonVariantSiteOnly_FirstSitesElementEmpty) {
    //one match crosses allele 1, and the other in nonvar
    const std::string prg_raw = "catttacatt5c6t5aaagcaacagaac";
    const uint64_t max_alphabet_num = 6;
    const std::string read = "acagaac";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    const Sites expected_sites = {
            {},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchOneNonVariantSiteOnly_DeleteFirstIntervalFalse) {
    //one match crosses allele 1, and the other in nonvar
    const std::string prg_raw = "catttacatt5c6t5aaagcaacagaac";
    const uint64_t max_alphabet_num = 6;
    const std::string read = "acagaac";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_FALSE(delete_first_interval);
}


TEST_F(BidirSearchBackward, MatchToMultipleNonVariantSitesOnly_SingleEmptySitesElement) {
    const std::string prg_raw = "catacagaacttacatt5g6t5aactagagagcaacagaactcacagaactc7cga8cgc8t";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "acagaac";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    const SA_Intervals expected_sa_intervals = {
            {6, 9},
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected_sites = {
            {},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchVariantSiteAndNonVariantSite) {
    //one match crosses allele 2, and the other in nonvar
    const std::string prg_raw = "catttacaca"
            "5g6t5"
            "aactagagagcaacataactctct";
    const uint64_t max_alphabet_num = 6;
    const std::string read = "acataac";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_FALSE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {5, 6},
            {6, 7},
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected_sites = {
            {},
            {VariantSite(5, {2})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchTwoLongVariantSites) {
    //read aligns from middle of  allele 3 of site 5 and allele 1 of site 7
    const std::string prg_raw = "acgacacat"
            "5gatag6tagga6gctcg6gctct5"
            "gctcgatgactagatagatag"
            "7cga8cgc8tga8tgc7"
            "ggcaacatctacga";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "gctcggctcgatgactagatagatagcgaggcaac";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 1, 1, 1, 1, 1,
            0, 2, 2, 2, 2, 2,
            0, 3, 3, 3, 3, 3,
            0, 4, 4, 4, 4, 4, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 1, 1, 1,
            0, 2, 2, 2,
            0, 3, 3, 3,
            0, 4, 4, 4, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_TRUE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {53, 54}
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected_sites = {
            {VariantSite(7, {1}), VariantSite(5, {})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, ReadStartsInFirstAllele_AlleleMissingFromSitesAlleleVector) {
    // Read aligns from middle of allele 3 of site 5 and allele 1 of site 7
    const std::string prg_raw = "acga"
            "5gctct6tt5"
            "gatat";
    const uint64_t max_alphabet_num = 6;
    const std::string read = "ctctgata";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 1, 1, 1, 1, 1,
            0, 2, 2, 0,
            0, 0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    const Sites expected_sites = {
            {VariantSite(5, {})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, ReadStartsInSecondAllele_AlleleMissingFromSitesAlleleVector) {
    // Read aligns from middle of allele 3 of site 5 and allele 1 of site 7
    const std::string prg_raw = "acga"
            "5tt6gctct5"
            "gatat";
    const uint64_t max_alphabet_num = 6;
    const std::string read = "ctctgata";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 1, 1,
            0, 2, 2, 2, 2, 2, 0,
            0, 0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    const Sites expected_sites = {
            {VariantSite(5, {})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, ReadEndsInSecondAllele_AlleleNumIncludedInSitesAlleleVector) {
    // Read aligns from middle of allele 3 of site 5 and allele 1 of site 7
    const std::string prg_raw = "acgc"
            "5tt6agata5"
            "tatag";
    const uint64_t max_alphabet_num = 6;
    const std::string read = "cgcagat";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 1, 1,
            0, 2, 2, 2, 2, 2, 0,
            0, 0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    const Sites expected_sites = {
            {VariantSite(5, {2})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchTwoVariantSites_FirstMatchVariantSiteHasEmptyAlleleVector) {
    // Read aligns from middle of allele 3 of site 5 and allele 1 of site 7
    const std::string prg_raw = "acgacacat"
            "5gatag6tagga6gctcg6gctct5"
            "gctcgatgactagatagatag"
            "7cga8cgc8tga8tgc7"
            "ggcaacatctacga";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "gctcggctcgatgactagatagatagcgaggcaac";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 0, 2, 2, 2, 2, 2, 0,
            3, 3, 3, 3, 3, 0, 4, 4, 4, 4, 4,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 1, 1, 0, 2, 2, 2, 0, 3, 3,
            3, 0, 4, 4, 4, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    const Sites expected_sites = {
            {VariantSite(7, {1}), VariantSite(5, {})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchWithinAlleleAndNonVariantSiteNoBoundaryCross_SitesVariantEmptyElement) {
    //read aligns in allele 2 of site 5, and in non-var region
    const std::string prg_raw = "gacatagacacacagt"
            "5gtcgcctcgtcggctttgagt6gtcgctgctccacacagagact5"
            "ggtgctagac"
            "7c8a7"
            "tcagctgctccacacagaga";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "ctgctccacacagaga";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0, 0,
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_FALSE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {45, 47}
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected_sites = {
            {}
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchWithinAlleleNoCrossingBoundary_SitesVariantEmptyElement) {
    //read aligns in allele 2 of site 5, and in non-var region
    const std::string prg_raw = "gacatagacacacagt"
            "5gtcgcctcgtcggctttgagt6gtcgctgctccacacagagact5"
            "ggtgctagac"
            "7c8a7"
            "tcag";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "ctgctccacacagaga";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0,
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_FALSE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {35, 36}
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected_sites = {
            {}
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchLongSiteRepeatedSnpOnSiteEdge) {
    //read aligns across sites 5 and 7, allele 1 in both cases
    const std::string prg_raw = "gacatagacacacagt"
            "5gtcgcctcgtcggctttgagt6gtcgctgctccacacagagact5"
            "ggtgctagac"
            "7c8a7"
            "ccagctgctccacacagaga";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "tagacacacagtgtcgcctcgtcggctttgagtggtgctagacccca";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0,
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_TRUE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {75, 76}
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected_sites = {
            {VariantSite(7, {1}), VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, MatchOverMultipleSites) {
    //read aligns over allele 1 of site 5, the nonVariantregion and allele 3 of site 7
    const std::string prg_raw = "acgacacat"
            "5gatag6tagga6gctcg6gctct5"
            "gctcgtgataatgactagatagatag"
            "7cga8cgc8tga8tgc7"
            "taggcaacatctacga";
    const uint64_t max_alphabet_num = 8;
    const std::string read = "tgata";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0,
            0, 1, 1, 1, 1, 1,
            0, 2, 2, 2, 2, 2,
            0, 3, 3, 3, 3, 3,
            0, 4, 4, 4, 4, 4, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0,
            0, 1, 1, 1,
            0, 2, 2, 2,
            0, 3, 3, 3,
            0, 4, 4, 4, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0,
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    //each element on the list corresponds to a SA interval
    //these elements are vectors of pairs (pair=(site, list of alleles))
    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_FALSE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {79, 80},
            {80, 81},
            {82, 83},
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    // note this unit test allows for an implementation limitation
    // of gramtools right now - unless a read crosses an odd number, it is not stored in sites()
    // should really notice the read has overlapped allele 3 of site 7, but it does not.
    const Sites expected_sites = {
            // first SA interval will be the match in the nonVariant region.
            // so we should get a vector of length zero, as it crosses no sites.
            {},

            // move to next SA interval - next element of list (sites)
            // this will be the overlap with site 7
            {VariantSite(7, {})},

            //next SA interval - overlap with site 5
            {VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected_sites);
}


TEST_F(BidirSearchBackward, SingleMatchOverManySites) {
    //overlaps site5-allele1, site7-allele2, site9-allele1, site11-allele1,  site13-allele2, site15-allele2
    const std::string prg_raw = "agggccta"
            "5c6t5"
            "acatgatc"
            "7a8g7"
            "tgatca"
            "9c10a9"
            "cata"
            "11g12t11"
            "aggtcgct"
            "13c14g13"
            "ggtc"
            "15atc16cat15"
            "ttcg";
    const uint64_t max_alphabet_num = 16;
    const std::string read = "cctacacatgatcgtgatcaccatagaggtcgctgggtccat";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 1,
            0, 2, 0,
            0, 0, 0, 0,
            0, 1, 1, 1,
            0, 2, 2, 2, 0,
            0, 0, 0, 0,
            0
    };

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval,
                     encoded_read.begin(), encoded_read.end(),
                     allele_mask, max_alphabet_num, kmer_index_generated,
                     rank_all, fm_index);

    EXPECT_TRUE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {
            {19, 20}
    };
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    // checking overlaps:
    // site5-allele1, site7-allele2, site9-allele1,
    // site11-allele1,  site13-allele2, site15-allele2
    const Sites expected_sites = {
            {
                    VariantSite(15, {2}),
                    VariantSite(13, {2}),
                    VariantSite(11, {1}),
                    VariantSite(9, {1}),
                    VariantSite(7, {2}),
                    VariantSite(5, {1}),
            }
    };
    EXPECT_EQ(sites, expected_sites);
}
