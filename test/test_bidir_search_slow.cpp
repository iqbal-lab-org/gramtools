#include <sdsl/suffix_arrays.hpp>
#include "gtest/gtest.h"

#include "process_prg.hpp"
#include "bwt_search.h"



void perform_test(const std::string &);

std::vector<string> generate_all_substrings(std::string);


TEST(SlowBackwardSearchTest, NoVariantsTest2) {
    std::string test_file2 = "./test_cases/13a.txt";
    perform_test(test_file2);
}


TEST(SlowBackwardSearchTest, NoVariantsABCABCTest3) {
    std::string test_file2 = "./test_cases/abc_abc_abc.txt";
    perform_test(test_file2);
}


TEST(SlowBackwardSearchTest, NoVariantsACTG4) {
    std::string test_file2 = "./test_cases/actg.txt";
    perform_test(test_file2);
}


TEST(SlowBackwardSearchTest, NoVariants_MSP34_200bp_Test5) {
    std::string test_file2 = "./test_cases/MSP3.4_200_bases.txt";
    perform_test(test_file2);
}


void perform_test(const std::string &test_fpath) {
    std::string prg;
    std::vector<string> substrings;
    std::vector<int> mask_a;

    //generate all substrings of PRG, use them all as queries
    ifstream ff(test_fpath);
    ff >> prg;
    substrings = generate_all_substrings(prg);

    //dummy mask
    uint32_t a;
    mask_a.clear();
    for (a = 0; a < prg.length(); a++) {
        mask_a.push_back(0);
    }

    FM_Index fm_index = construct_fm_index(test_fpath,
                                           "int_alphabet_file",
                                           "memory_log_file",
                                           "csa_file", true);

    std::list<std::pair<uint64_t, uint64_t>> sa_intervals, sa_intervals_rev;
    std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;

    std::vector<uint8_t> p_tmp;
    std::string q_tmp;

    for (std::vector<std::string>::iterator it = substrings.begin(); it < substrings.end(); ++it) {
        q_tmp = *it;

        bool first_del = false;
        bool precalc = false;
        int occ_expt = 0;
        auto pos = prg.find(q_tmp, 0);

        while (pos != std::string::npos) {
            occ_expt++;
            pos = prg.find(q_tmp, pos + 1);
        }

        for (uint16_t i = 0; i < q_tmp.length(); i++) {
            if (q_tmp[i] == 'A' or q_tmp[i] == 'a') p_tmp.push_back(1);
            if (q_tmp[i] == 'C' or q_tmp[i] == 'c') p_tmp.push_back(2);
            if (q_tmp[i] == 'G' or q_tmp[i] == 'g') p_tmp.push_back(3);
            if (q_tmp[i] == 'T' or q_tmp[i] == 't') p_tmp.push_back(4);
        }

        bidir_search_bwd(fm_index, 0, fm_index.size(), 0, fm_index.size(), p_tmp.begin(),
                         p_tmp.end(),
                         sa_intervals, sa_intervals_rev, sites, mask_a, 5,
                         first_del, precalc);

        uint64_t no_occ = (*sa_intervals.begin()).second - (*sa_intervals.begin()).first;
        EXPECT_FALSE(first_del);
        EXPECT_EQ(1, sa_intervals.size());
        EXPECT_EQ(no_occ, occ_expt);

        sa_intervals.clear();
        sa_intervals_rev.clear();
        sites.clear();

        p_tmp.clear();
    }
}


std::vector<std::string> generate_all_substrings(std::string q) {
    std::vector<std::string> substrings;
    auto c = 0, i = 0;

    auto n = q.size();
    for (c = 0; c < n; c++) {
        for (i = 1; i <= n - c; i++) {
            substrings.push_back(q.substr(c, c + i));
        }
    }
    return substrings;
}
