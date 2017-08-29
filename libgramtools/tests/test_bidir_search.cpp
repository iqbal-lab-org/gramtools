#include <iostream>
#include <vector>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <sdsl/suffix_arrays.hpp>

#include "gtest/gtest.h"

#include "common.hpp"
#include "map.hpp"
#include "bidir_search_bwd.hpp"


void set_up_state(FM_Index &fm_index, DNA_Rank &rank_all, std::vector<int> &allele_mask,
                  const std::string &prg_fpath, const std::string &mask_fpath) {
    std::ifstream g(mask_fpath);
    int a;
    while (g >> a)
        allele_mask.push_back(a);

    fm_index = get_fm_index(true, "csa_file", "int_alphabet_file", prg_fpath, "memory_log_file");
    rank_all = calculate_ranks(fm_index);
}


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


class BackwardSearchTestFileFree : public ::testing::Test {

protected:
    std::string prg_fpath;

    virtual void SetUp() {
        boost::uuids::uuid uuid = boost::uuids::random_generator()();
        const auto uuid_str = boost::lexical_cast<std::string>(uuid);
        prg_fpath = "./prg_" + uuid_str;
    }

    FM_Index fm_index_from_raw_prg(const std::string &prg_raw) {
        std::vector<uint64_t> prg = encode_prg(prg_raw);
        dump_encoded_prg(prg, prg_fpath);
        FM_Index fm_index;
        sdsl::construct(fm_index, prg_fpath, 8);
        return fm_index;
    }

    virtual void TearDown() {
        std::remove(prg_fpath.c_str());
    }

};


TEST_F(BackwardSearchTestFileFree, NoVariants1) {
    const auto prg_raw = "a";
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
                     allele_mask, 4, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = (*sa_intervals.begin()).second - (*sa_intervals.begin()).first;
    EXPECT_TRUE(!delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 1);
    EXPECT_EQ(no_occ, 1);

    const Sites expected = {Site()};
    EXPECT_EQ(expected, sites);
}


TEST_F(BackwardSearchTestFileFree, OneSNP) {
    // aligns across SNP allele 1 (and both flanks)
    const auto prg_raw = "catttacaca5g6t5aactagagagca";
    const auto read = "ttacacagaactagagag";
    const std::vector<int> allele_mask = {0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 1, 0, 2, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0};

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_read(read);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 6, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = (*sa_intervals.begin()).second - (*sa_intervals.begin()).first;
    EXPECT_TRUE(delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 1);
    EXPECT_EQ(no_occ, 1);

    const Sites expected = {
            {VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
}


TEST(BackwardSearchTest, TwoSNPs) {
    cleanup_files();

    //prg = catttacaca5g6t5aactag7a8g7agcagggt
    const auto prg_fpath = "./test_cases/two_snps.txt";
    const auto mask_fpath = "./test_cases/two_snps_mask_a.txt";

    //aligns across both SNPs, both allele 1
    const auto read = "ttacacagaactagaagcag";
    const auto encoded_read = encode_read(read);

    FM_Index fm_index;
    DNA_Rank rank_all;
    std::vector<int> allele_mask;
    set_up_state(fm_index, rank_all, allele_mask, prg_fpath, mask_fpath);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 8, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = (*sa_intervals.begin()).second - (*sa_intervals.begin()).first;
    EXPECT_TRUE(delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 1);
    EXPECT_EQ(no_occ, 1);

    const Sites expected = {
            {VariantSite(7, {1}), VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
    cleanup_files();
}


TEST(BackwardSearchTest, Two_matches_one_variable_one_nonvariable_region) {
    cleanup_files();

    //prg=catttacaca5g6t5aactagagagcaacagaactctct
    const auto prg_fpath = "./test_cases/two_matches_var_nonvar.txt";
    const auto mask_fpath = "./test_cases/two_matches_var_nonvar_mask_a.txt";

    //one match crosses allele 1, and the other in nonvar
    const auto read = "acagaac";
    const auto encoded_read = encode_read(read);

    std::ifstream g(mask_fpath);
    std::vector<int> allele_mask;
    int a;
    while (g >> a)
        allele_mask.push_back(a);

    const FM_Index fm_index = get_fm_index(true, "csa_file", "int_alphabet_file", prg_fpath, "memory_log_file");
    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 6, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = 0;
    for (const auto &sa_interval: sa_intervals)
        no_occ += sa_interval.second - sa_interval.first;

    EXPECT_TRUE(!delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 2);
    EXPECT_EQ(no_occ, 2);

    const Sites expected = {
            {},
            {VariantSite(5, {1})},
    };
    EXPECT_EQ(sites, expected);
    cleanup_files();
}


TEST(BackwardSearchTest, Two_matches_one_variable_second_allele_one_nonvariable_region) {
    cleanup_files();

    //prg=catttacaca5g6t5aactagagagcaacataactctct
    const auto prg_fpath = "./test_cases/two_matches_var_other_allele_nonvar.txt";
    const auto mask_fpath = "./test_cases/two_matches_var_nonvar_mask_a.txt";

    //one match crosses allele 2, and the other in nonvar
    const auto read = "acataac";
    const auto encoded_read = encode_read(read);

    std::ifstream g(mask_fpath);
    std::vector<int> allele_mask;
    int a;
    while (g >> a)
        allele_mask.push_back(a);

    const FM_Index fm_index = get_fm_index(true, "csa_file", "int_alphabet_file", prg_fpath, "memory_log_file");
    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 6, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = 0;
    for (const auto &sa_interval: sa_intervals)
        no_occ += sa_interval.second - sa_interval.first;

    EXPECT_TRUE(!delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 2);
    EXPECT_EQ(no_occ, 2);

    const Sites expected = {
            {},
            {VariantSite(5, {2})},
    };
    EXPECT_EQ(sites, expected);
    cleanup_files();
}


TEST(BackwardSearchTest, Two_long_sites) {
    cleanup_files();
    //prg = acgacacat5gatag6tagga6gctcg6gctct5gctcgatgactagatagatag7cga8cgc8tga8tgc7ggcaacatctacga
    const auto prg_fpath = "./test_cases/two_long_sites.txt";
    const auto mask_fpath = "./test_cases/two_long_sites_mask_a.txt";

    //read aligns from middle of  allele 3 of site 5 and allele 1 of site 7
    const auto read = "gctcggctcgatgactagatagatagcgaggcaac";
    const auto encoded_read = encode_read(read);


    std::ifstream g(mask_fpath);
    std::vector<int> allele_mask;
    int a;
    while (g >> a)
        allele_mask.push_back(a);

    const FM_Index fm_index = get_fm_index(true, "csa_file", "int_alphabet_file", prg_fpath, "memory_log_file");
    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 8, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = 0;
    for (const auto &sa_interval: sa_intervals)
        no_occ += sa_interval.second - sa_interval.first;

    EXPECT_TRUE(delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 1);
    EXPECT_EQ(no_occ, 1);

    const Sites expected = {
            {VariantSite(7, {1}), VariantSite(5, {})},
    };
    EXPECT_EQ(sites, expected);
    cleanup_files();
}


TEST(BackwardSearchTest, Match_within_long_site_match_outside) {
    cleanup_files();

    //prg=gacatagacacacagt5gtcgcctcgtcggctttgagt6gtcgctgctccacacagagact5ggtgctagac7c8a7tcagctgctccacacagaga
    const auto prg_fpath = "./test_cases/match_within_long_site.txt";
    const auto mask_fpath = "./test_cases/match_within_long_site_mask_a.txt";

    //read aligns in allele 2 of site 5, and in non-var region
    const auto read = "ctgctccacacagaga";
    const auto encoded_read = encode_read(read);

    std::ifstream g(mask_fpath);
    std::vector<int> allele_mask;
    int a;
    while (g >> a)
        allele_mask.push_back(a);

    const FM_Index fm_index = get_fm_index(true, "csa_file", "int_alphabet_file", prg_fpath, "memory_log_file");
    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 8, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = 0;
    for (const auto &sa_interval: sa_intervals)
        no_occ += sa_interval.second - sa_interval.first;

    EXPECT_TRUE(!delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 1);
    EXPECT_EQ(no_occ, 2);

    const Sites expected = {Site()};
    EXPECT_EQ(expected, sites);
    cleanup_files();
}


TEST(BackwardSearchTest, Long_site_and_repeated_snp_on_edge_of_site) {
    cleanup_files();

    //prg = gacatagacacacagt5gtcgcctcgtcggctttgagt6gtcgctgctccacacagagact5ggtgctagac7c8a7ccagctgctccacacagaga
    const auto prg_fpath = "./test_cases/repeated_snp_on_both_edges.txt";
    const auto mask_fpath = "./test_cases/match_within_long_site_mask_a.txt";

    //read aligns across sites 5 and 7, allele 1 in both cases
    const auto read = "tagacacacagtgtcgcctcgtcggctttgagtggtgctagacccca";
    const auto encoded_read = encode_read(read);

    std::ifstream g(mask_fpath);
    std::vector<int> allele_mask;
    int a;
    while (g >> a)
        allele_mask.push_back(a);

    const FM_Index fm_index = get_fm_index(true, "csa_file", "int_alphabet_file", prg_fpath, "memory_log_file");
    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 8, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = 0;
    for (const auto &sa_interval: sa_intervals)
        no_occ += sa_interval.second - sa_interval.first;

    EXPECT_TRUE(delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 1);
    EXPECT_EQ(no_occ, 1);

    const Sites expected = {
            {VariantSite(7, {1}), VariantSite(5, {1})},
    };
    EXPECT_EQ(expected, sites);
    cleanup_files();
}


TEST(BackwardSearchTest, Multiple_matches_over_multiple_sites) {
    cleanup_files();

    //prg=acgacacat5gatag6tagga6gctcg6gctct5gctcgtgataatgactagatagatag7cga8cgc8tga8tgc7taggcaacatctacga
    const auto prg_fpath = "./test_cases/multiple_matches_multiple_sites.txt";
    const auto mask_fpath = "./test_cases/multiple_matches_multiple_sites_mask_a.txt";

    //read aligns over allele 1 of site 5, the nonvariableregion and allele 3 of site 7
    const auto read = "tgata";
    const auto encoded_read = encode_read(read);

    std::ifstream g(mask_fpath);
    std::vector<int> allele_mask;
    int a;
    while (g >> a)
        allele_mask.push_back(a);

    const FM_Index fm_index = get_fm_index(true, "csa_file", "int_alphabet_file", prg_fpath, "memory_log_file");
    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    //each element on the list corresponds to a SA interval
    //these elements are vectors of pairs (pair=(site, list of alleles))
    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 8, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = 0;
    for (const auto &sa_interval: sa_intervals)
        no_occ += sa_interval.second - sa_interval.first;

    EXPECT_FALSE(delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 3);
    EXPECT_EQ(no_occ, 3);

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
    cleanup_files();
}


TEST(BackwardSearchTest, One_match_many_sites) {
    cleanup_files();

    //prg=agggccta5c6t5acatgatc7a8g7tgatca9c10a9cata11g12t11aggtcgct13c14g13ggtc15atc16cat15ttcg
    const auto prg_fpath = "./test_cases/One_match_many_sites.txt";
    const auto mask_fpath = "./test_cases/One_match_many_sites_mask_a.txt";

    //overlaps site5-allele1, site7-allele2, site9-allele1, site11-allele1,  site13-allele2, site15-allele2
    const auto read = "cctacacatgatcgtgatcaccatagaggtcgctgggtccat";
    const auto encoded_read = encode_read(read);

    std::ifstream g(mask_fpath);
    std::vector<int> allele_mask;
    int a;
    while (g >> a)
        allele_mask.push_back(a);

    const FM_Index fm_index = get_fm_index(true, "csa_file", "int_alphabet_file", prg_fpath, "memory_log_file");
    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    Sites sites = {Site()};

    bool delete_first_interval = false;
    const bool kmer_index_generated = false;

    bidir_search_bwd(sa_intervals, sites, delete_first_interval, encoded_read.begin(),
                     encoded_read.end(), allele_mask, 16, kmer_index_generated, rank_all, fm_index);

    uint64_t no_occ = 0;
    for (const auto &sa_interval: sa_intervals)
        no_occ += sa_interval.second - sa_interval.first;

    EXPECT_TRUE(delete_first_interval);
    EXPECT_EQ(sa_intervals.size(), 1);
    EXPECT_EQ(no_occ, 1);

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
    cleanup_files();
}
