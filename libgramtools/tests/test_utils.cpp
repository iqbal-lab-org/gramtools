#include "test_utils.hpp"


PRG_Info generate_prg_info(const std::string &prg_raw) {
    Parameters parameters;
    parameters.encoded_prg_fpath = "@encoded_prg_file_name";
    parameters.fm_index_fpath = "@fm_index";

    auto encoded_prg = encode_prg(prg_raw);
    sdsl::store_to_file(encoded_prg, parameters.encoded_prg_fpath);

    PRG_Info prg_info;
    prg_info.fm_index = generate_fm_index(parameters);
    prg_info.sites_mask = generate_sites_mask(prg_raw);
    prg_info.allele_mask = generate_allele_mask(prg_raw);
    prg_info.max_alphabet_num = max_alphabet_num(prg_raw);

    // prg_info.dna_rank = calculate_ranks(prg_info.fm_index);
    return prg_info;
}
