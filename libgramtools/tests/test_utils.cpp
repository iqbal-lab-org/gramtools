#include "gtest/gtest.h"

#include "masks.hpp"
#include "kmer_index.hpp"
#include "test_utils.hpp"


std::pair<uint64_t, std::string::const_iterator> get_marker(const std::string::const_iterator start_it,
                                                            const std::string::const_iterator end_it) {
    std::string::const_iterator it = start_it;
    std::string digits;
    while (isdigit(*it) and it != end_it) {
        digits += *it;
        it++;
    }
    it--;
    uint64_t marker = (uint64_t) std::atoi(digits.c_str());
    return std::make_pair(marker, it);
}


uint64_t max_alphabet_num(const std::string &prg_raw) {
    uint64_t max_alphabet_num = 1;

    auto it = prg_raw.begin();
    auto end_it = prg_raw.end();
    while (it != prg_raw.end()) {
        const auto not_marker = !isdigit(*it);
        if (not_marker) {
            uint64_t encoded_base = (uint64_t) encode_dna_base(*it);
            if (encoded_base > max_alphabet_num)
                max_alphabet_num = encoded_base;
        } else {
            uint64_t marker = 0;
            std::tie(marker, it) = get_marker(it, end_it);
            if (marker > max_alphabet_num)
                max_alphabet_num = marker;
        }
        it++;
    }
    return max_alphabet_num;
}


std::vector<Marker> generate_sites_mask(const std::string &prg_raw) {
    std::vector<Marker> sites_mask;
    uint64_t current_site_edge_marker = 0;

    auto it = prg_raw.begin();
    auto end_it = prg_raw.end();
    while (it != prg_raw.end()) {
        const auto not_marker = not isdigit(*it);
        if (not_marker) {
            sites_mask.push_back(current_site_edge_marker);
            it++;
            continue;
        }

        uint64_t marker = 0;
        std::tie(marker, it) = get_marker(it, end_it);
        it++;

        sites_mask.push_back(0);

        bool site_edge_marker = marker % 2 != 0;
        if (site_edge_marker) {
            const auto at_site_start = current_site_edge_marker == 0;
            if (at_site_start) {
                current_site_edge_marker = marker;
            } else {
                current_site_edge_marker = 0;
            }
            continue;
        }
    }
    return sites_mask;
}


PRG_Info generate_prg_info(const std::string &prg_raw) {
    Parameters parameters;
    parameters.encoded_prg_fpath = "@encoded_prg_file_name";
    parameters.fm_index_fpath = "@fm_index";

    auto encoded_prg = encode_prg(prg_raw);
    sdsl::store_to_file(encoded_prg, parameters.encoded_prg_fpath);

    PRG_Info prg_info;
    prg_info.fm_index = generate_fm_index(parameters);
    prg_info.encoded_prg = encoded_prg;
    prg_info.sites_mask = generate_sites_mask(prg_raw);
    prg_info.allele_mask = generate_allele_mask(encoded_prg);
    prg_info.max_alphabet_num = max_alphabet_num(prg_raw);

    // prg_info.dna_rank = calculate_ranks(prg_info.fm_index);
    return prg_info;
}


TEST(GenerateSitesMask, SingleVariantSiteTwoAlleles_CorrectSitesMask) {
    const std::string prg_raw = "a5g6t5c";
    auto result = generate_sites_mask(prg_raw);
    std::vector<Marker> expected = {
            0, 0, 5, 0, 5, 0, 0
    };
    EXPECT_EQ(result, expected);
}


TEST(GenerateSitesMask, TwoVariantSites_CorrectSitesMask) {
    const std::string prg_raw = "a5g6t5cc7g8tt8aa7";
    auto result = generate_sites_mask(prg_raw);
    std::vector<Marker> expected = {
            0,
            0, 5, 0, 5, 0,
            0, 0,
            0, 7, 0, 7, 7, 0, 7, 7, 0
    };
    EXPECT_EQ(result, expected);
}


TEST(MaxAlphabetNum, PrgWithVariantSite_LargestSiteMarkerAsMaxAlphabet) {
    const std::string prg_raw = "a13g14t13tt";
    auto result = max_alphabet_num(prg_raw);
    uint64_t expected = 14;
    EXPECT_EQ(result, expected);
}

TEST(MaxAlphabetNum, SingleCharPrg_CorrectBaseEncodingAsMaxAlphabet) {
    const std::string prg_raw = "c";
    auto result = max_alphabet_num(prg_raw);
    uint64_t expected = 2;
    EXPECT_EQ(result, expected);
}
