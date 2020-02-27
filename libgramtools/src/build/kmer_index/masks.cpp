#include <cstdint>
#include <vector>
#include <string>

#include "build/kmer_index/masks.hpp"


using namespace gram;

sdsl::int_vector<> gram::generate_allele_mask(const marker_vec &encoded_prg) {
    sdsl::int_vector<> allele_mask(encoded_prg.size(), 0, 32);
    uint32_t current_allele_id = 0;
    uint64_t last_allele_position = 0;

    for (uint64_t i = 0; i < encoded_prg.size(); ++i) {
        const auto &prg_char = encoded_prg[i];
        bool markable = prg_char <= 4;

        if (markable){
            if (current_allele_id > 0) allele_mask[i] = current_allele_id;
            continue;
        }

        // Implicitly, from here, we are dealing with a marker (prg_char > 4).
        bool is_site_marker = prg_char %2 == 1;
        if (is_site_marker){
            current_allele_id = 1;
            if (last_allele_position > 0){ // Is this not the first site we are seeing?
                // Reset to 0 those positions in between the last allele marker and the new variant site.
                for (int pos = last_allele_position + 1; pos < i; ++pos) allele_mask[pos] = 0;
            }
            continue;
        }
        else {
            ++current_allele_id;
            last_allele_position = i;
        }
    }
    // Reset to 0 those positions in between the last allele marker and the end of the PRG string.
    for (int pos = last_allele_position + 1; pos < encoded_prg.size(); ++pos) allele_mask[pos] = 0;
    sdsl::util::bit_compress(allele_mask);
    return allele_mask;
}


sdsl::int_vector<> gram::load_allele_mask(CommonParameters const &parameters) {
    sdsl::int_vector<> allele_mask;
    sdsl::load_from_file(allele_mask, parameters.allele_mask_fpath);
    return allele_mask;
}


sdsl::int_vector<> gram::generate_sites_mask(const marker_vec &encoded_prg) {
    sdsl::int_vector<> sites_mask(encoded_prg.size(), 0, 32);
    Marker current_site_marker = 0;
    uint64_t last_allele_position = 0;

    for (uint64_t i = 0; i < encoded_prg.size(); ++i) {
        const auto &prg_char = encoded_prg[i];
        bool markable = prg_char <= 4;

        if (markable){
            if (current_site_marker > 0) sites_mask[i] = current_site_marker;
            continue;
        }

        // Implicitly, from here, we are dealing with a marker (prg_char > 4).
        bool is_site_marker = prg_char %2 == 1;
        if (is_site_marker){
            current_site_marker = prg_char;
            if (last_allele_position > 0){ // Is this not the first site we are seeing?
                // Reset to 0 those positions in between the last allele marker and the new variant site.
                for (int pos = last_allele_position + 1; pos < i; ++pos) sites_mask[pos] = 0;
            }
            continue;
        }
        else {last_allele_position = i;}
    }
    // Reset to 0 those positions in between the last allele marker and the end of the PRG string.
    for (int pos = last_allele_position + 1; pos < encoded_prg.size(); ++pos) sites_mask[pos] = 0;

    sdsl::util::bit_compress(sites_mask);
    return sites_mask;
}

sdsl::int_vector<> gram::load_sites_mask(CommonParameters const &parameters) {
    sdsl::int_vector<> sites_mask;
    sdsl::load_from_file(sites_mask, parameters.sites_mask_fpath);
    return sites_mask;
}

sdsl::bit_vector gram::generate_prg_markers_mask(const marker_vec &encoded_prg) {
    sdsl::bit_vector variants_markers_mask(encoded_prg.size(), 0);
    for (uint64_t i = 0; i < encoded_prg.size(); i++)
        variants_markers_mask[i] = encoded_prg[i] > 4;
    return variants_markers_mask;
}
