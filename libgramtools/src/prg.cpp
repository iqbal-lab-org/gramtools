#include "kmers.hpp"
#include "prg.hpp"


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


std::vector<int> generate_allele_mask(const std::string &prg_raw) {
    std::vector<int> allele_mask;
    uint64_t current_site_edge_marker = 0;
    int current_allele_number = 0;

    auto it = prg_raw.begin();
    auto end_it = prg_raw.end();
    while (it != prg_raw.end()) {
        const auto not_marker = !isdigit(*it);
        if (not_marker) {
            allele_mask.push_back(current_allele_number);
            it++;
            continue;
        }

        uint64_t marker = 0;
        std::tie(marker, it) = get_marker(it, end_it);
        it++;

        allele_mask.push_back(0);

        bool site_edge_marker = marker % 2 != 0;
        if (site_edge_marker) {
            const auto at_site_start = current_site_edge_marker == 0;
            if (at_site_start) {
                current_site_edge_marker = marker;
                current_allele_number = 1;
            } else {
                current_site_edge_marker = 0;
                current_allele_number = 0;
            }
            continue;
        } else {
            current_allele_number++;
        }
    }
    return allele_mask;
}