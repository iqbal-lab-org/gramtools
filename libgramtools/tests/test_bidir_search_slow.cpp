#include <sdsl/suffix_arrays.hpp>
#include "gtest/gtest.h"

#include "map.hpp"
#include "bidir_search_bwd.hpp"


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


void perform_test(const std::string &test_fpath) {
    std::string prg;
    std::vector<std::string> substrings;
    std::vector<int> allele_mask;

    //generate all substrings of PRG, use them all as queries
    std::ifstream ff(test_fpath);
    ff >> prg;
    substrings = generate_all_substrings(prg);

    //dummy mask
    uint64_t a;
    allele_mask.clear();
    for (a = 0; a < prg.length(); a++) {
        allele_mask.push_back(0);
    }

    const FM_Index fm_index = construct_fm_index(true, "csa_file", "int_alphabet_file", test_fpath, "memory_log_file");

    SA_Intervals sa_intervals;
    Sites sites;

    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    std::vector<uint8_t> p_tmp;
    std::string q_tmp;

    for (auto it = substrings.begin(); it < substrings.end(); ++it) {
        q_tmp = *it;

        bool delete_first = false;
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

        if (sa_intervals.empty()) {
            sa_intervals.emplace_back(std::make_pair(0, fm_index.size()));

            Site empty_pair_vector;
            sites.push_back(empty_pair_vector);
        }

        bidir_search_bwd(sa_intervals, sites, delete_first, p_tmp.begin(),
                         p_tmp.end(), allele_mask, 5, precalc, rank_all, fm_index);

        uint64_t no_occ = (*sa_intervals.begin()).second - (*sa_intervals.begin()).first;
        EXPECT_FALSE(delete_first);
        EXPECT_EQ(1, sa_intervals.size());
        EXPECT_EQ(no_occ, occ_expt);

        sa_intervals.clear();
        sites.clear();

        p_tmp.clear();
    }
}


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
