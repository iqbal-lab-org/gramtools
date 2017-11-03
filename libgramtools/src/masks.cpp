#include <cstdint>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "parameters.hpp"
#include "masks.hpp"


sdsl::bit_vector generate_markers_mask(const sdsl::int_vector<> &encoded_prg) {
    sdsl::bit_vector variants_markers_mask(encoded_prg.size(), 0);
    for (unsigned int i = 0; i < encoded_prg.size(); i++)
        variants_markers_mask[i] = encoded_prg[i] > 4;
    return variants_markers_mask;
}


sdsl::int_vector<> load_allele_mask(const Parameters &parameters) {
    sdsl::int_vector<> allele_mask;
    sdsl::load_from_file(allele_mask, parameters.allele_mask_fpath);
    return allele_mask;
}


sdsl::int_vector<> generate_allele_mask(const sdsl::int_vector<> &encoded_prg) {
    sdsl::int_vector<> allele_mask(encoded_prg.size(), 0, 32);
    uint32_t current_allele_id = 1;
    bool within_variant_site = false;

    for (uint64_t i = 0; i < encoded_prg.size(); ++i) {
        const auto &prg_char = encoded_prg[i];
        auto at_varaint_site_boundary = prg_char > 4
                                        and prg_char % 2 != 0;
        auto entering_variant_site = at_varaint_site_boundary
                                     and not within_variant_site;
        if (entering_variant_site){
            within_variant_site = true;
            current_allele_id = 1;
            continue;
        }

        auto exiting_variant_site = at_varaint_site_boundary
                                     and within_variant_site;
        if (exiting_variant_site){
            within_variant_site = false;
            continue;
        }

        auto within_allele = prg_char <= 4
                             and within_variant_site;
        if (within_allele) {
            allele_mask[i] = current_allele_id;
            continue;
        }

        auto at_allele_marker = prg_char > 4 and prg_char % 2 == 0;
        if (at_allele_marker) {
            current_allele_id++;
            continue;
        }
    }
    return allele_mask;
}


MasksParser::MasksParser(const std::string &sites_fname) {
    std::ifstream sites_stream(sites_fname);
    parse_sites(sites_stream);
}


void MasksParser::parse_sites(std::istream &stream) {
    uint64_t max_sites_count = 0;
    uint64_t site_count;

    while (stream >> site_count) {
        MasksParser::sites.push_back(site_count);
        if (site_count > max_sites_count)
            max_sites_count = site_count;
    }

    // no_sites is last odd number in mask_sites, but alphabet size
    // is the even number corresponding to it
    MasksParser::max_alphabet_num = max_sites_count + 1;
}
