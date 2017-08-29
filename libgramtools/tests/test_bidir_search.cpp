#include <iostream>
#include <vector>
#include <fstream>

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


class BackwardSearchTest : public ::testing::Test {

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


TEST_F(BackwardSearchTest, NoVariants1) {
    const auto prg_raw = "a";
    const uint64_t max_alphabet_num = 4;
    const std::string read = "a";
    std::vector<int> allele_mask = {0};

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

    const SA_Intervals expected_sa_intervals = {{1, 2}};
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected = {Site()};
    EXPECT_EQ(expected, sites);
}


TEST_F(BackwardSearchTest, OneSNP) {
    // aligns across SNP allele 1 (and both flanks)
    const auto prg_raw = "catttacaca5g6t5aactagagagca";
    const uint64_t max_alphabet_num = 6;
    const auto read = "ttacacagaactagagag";
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

    const SA_Intervals expected_sa_intervals = {{22, 23}};
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected = {
            {VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BackwardSearchTest, TwoSNPs) {
    //aligns across both SNPs, both allele 1
    const auto prg_raw = "catttacaca5g6t5aactag7a8g7agcagggt";
    const uint64_t max_alphabet_num = 8;
    const auto read = "ttacacagaactagaagcag";
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

    const SA_Intervals expected_sa_intervals = {{26, 27}};
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected = {
            {VariantSite(7, {1}), VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BackwardSearchTest, Two_matches_one_variable_one_nonvariable_region) {
    //one match crosses allele 1, and the other in nonvar
    const auto prg_raw = "catttacaca5g6t5aactagagagcaacagaactctct";
    const uint64_t max_alphabet_num = 6;
    const auto read = "acagaac";
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

    const Sites expected = {
            {},
            {VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BackwardSearchTest, Two_matches_one_variable_second_allele_one_nonvariable_region) {
    //one match crosses allele 2, and the other in nonvar
    const auto prg_raw = "catttacaca5g6t5aactagagagcaacataactctct";
    const uint64_t max_alphabet_num = 6;
    const auto read = "acataac";
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

    const Sites expected = {
            {},
            {VariantSite(5, {2})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BackwardSearchTest, Two_long_sites) {
    //read aligns from middle of  allele 3 of site 5 and allele 1 of site 7
    const auto prg_raw = "acgacacat5gatag6tagga6gctcg6gctct5gctcgatgactagatagatag7cga8cgc8tga8tgc7ggcaacatctacga";
    const uint64_t max_alphabet_num = 8;
    const auto read = "gctcggctcgatgactagatagatagcgaggcaac";
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

    EXPECT_TRUE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {{53, 54}};
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected = {
            {VariantSite(7, {1}), VariantSite(5, {})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BackwardSearchTest, Match_within_long_site_match_outside) {
    //read aligns in allele 2 of site 5, and in non-var region
    const auto prg_raw = "gacatagacacacagt5gtcgcctcgtcggctttgagt6gtcgctgctccacacagagact5ggtgctagac7c8a7tcagctgctccacacagaga";
    const uint64_t max_alphabet_num = 8;
    const auto read = "ctgctccacacagaga";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 0, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 2, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
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

    const SA_Intervals expected_sa_intervals = {{45, 47}};
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected = {Site()};
    EXPECT_EQ(expected, sites);
}


TEST_F(BackwardSearchTest, Long_site_and_repeated_snp_on_edge_of_site) {
    //read aligns across sites 5 and 7, allele 1 in both cases
    const auto prg_raw = "gacatagacacacagt5gtcgcctcgtcggctttgagt6gtcgctgctccacacagagact5ggtgctagac7c8a7ccagctgctccacacagaga";
    const uint64_t max_alphabet_num = 8;
    const auto read = "tagacacacagtgtcgcctcgtcggctttgagtggtgctagacccca";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 0, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 2, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
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

    const SA_Intervals expected_sa_intervals = {{75, 76}};
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const Sites expected = {
            {VariantSite(7, {1}), VariantSite(5, {1})},
    };
    EXPECT_EQ(expected, sites);
}


TEST_F(BackwardSearchTest, Multiple_matches_over_multiple_sites) {
    //read aligns over allele 1 of site 5, the nonvariableregion and allele 3 of site 7
    const auto prg_raw = "acgacacat5gatag6tagga6gctcg6gctct5gctcgtgataatgactagatagatag7cga8cgc8tga8tgc7taggcaacatctacga";
    const uint64_t max_alphabet_num = 8;
    const auto read = "tgata";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 1, 1, 1, 1, 0, 2, 2,
            2, 2, 2, 0, 3, 3, 3, 3, 3,
            0, 4, 4, 4, 4, 4, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 1,
            1, 0, 2, 2, 2, 0, 3, 3, 3,
            0, 4, 4, 4, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
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
    const Sites expected = {
            // first SA interval will be the match in the nonvariable region.
            // so we should get a vector of length zero, as it crosses no sites.
            {},

            // move to next SA interval - next element of list (sites)
            // this will be the overlap with site 7
            {VariantSite(7, {})},

            //next SA interval - overlap with site 5
            {VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
}


TEST_F(BackwardSearchTest, One_match_many_sites) {
    //overlaps site5-allele1, site7-allele2, site9-allele1, site11-allele1,  site13-allele2, site15-allele2
    const auto prg_raw = "agggccta5c6t5acatgatc7a8g7tgatca9c10a9cata11g12t11aggtcgct13c14g13ggtc15atc16cat15ttcg";
    const uint64_t max_alphabet_num = 16;
    const auto read = "cctacacatgatcgtgatcaccatagaggtcgctgggtccat";
    const std::vector<int> allele_mask = {
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 2, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 2, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 2,
            0, 0, 0, 0, 0, 0, 1, 0, 2,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 2, 0, 0, 0, 0, 0,
            0, 1, 1, 1, 0, 2, 2, 2, 0,
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

    EXPECT_TRUE(delete_first_interval);

    const SA_Intervals expected_sa_intervals = {{19, 20}};
    EXPECT_EQ(sa_intervals, expected_sa_intervals);

    const auto &first_site = sites.front();

    // checking overlaps:
    // site5-allele1, site7-allele2, site9-allele1,
    // site11-allele1,  site13-allele2, site15-allele2
    const Site expected_site = {
            VariantSite(15, {2}),
            VariantSite(13, {2}),
            VariantSite(11, {1}),
            VariantSite(9, {1}),
            VariantSite(7, {2}),
            VariantSite(5, {1}),
    };
    EXPECT_EQ(first_site, expected_site);
}
