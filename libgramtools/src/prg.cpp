#include "prg.hpp"


std::vector<int> generate_allele_mask(const std::string &prg_raw) {
    std::vector<int> allele_mask;

    int current_site_edge_marker = 0;
    int current_allele_number = 0;
    for (const auto &c: prg_raw) {
        const auto not_marker = !isdigit(c);
        if (not_marker) {
            allele_mask.push_back(current_allele_number);
            continue;
        }

        allele_mask.push_back(0);
        int marker = atoi(&c);

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