#include <boost/filesystem.hpp>

#include "utils.hpp"
#include "dna_ranks.hpp"


namespace fs = boost::filesystem;


sdsl::bit_vector generate_base_bwt_mask(const Base &base,
                                        const FM_Index &fm_index) {
    sdsl::bit_vector mask(fm_index.bwt.size(), 0);
    for (uint64_t i = 0; i < fm_index.bwt.size(); i++)
        mask[i] = fm_index.bwt[i] == base;
    return mask;
}


std::string bwt_mask_fname(const std::string &base_char,
                           const Parameters &parameters) {
    auto handling_unit_tests = parameters.gram_dirpath[0] == '@';
    if (handling_unit_tests) {
        return parameters.gram_dirpath
               + "_"
               + base_char
               + "_base_bwt_mask";
    }
    fs::path dir(parameters.gram_dirpath);
    fs::path file(base_char + "_base_bwt_mask");
    fs::path full_path = dir / file;
    return full_path.string();
}


void generate_dna_bwt_masks(const FM_Index &fm_index,
                            const Parameters &parameters) {
    auto a_mask = generate_base_bwt_mask(1, fm_index);
    auto c_mask = generate_base_bwt_mask(2, fm_index);
    auto g_mask = generate_base_bwt_mask(3, fm_index);
    auto t_mask = generate_base_bwt_mask(4, fm_index);

    auto fpath = bwt_mask_fname("a", parameters);
    sdsl::store_to_file(a_mask, fpath);

    fpath = bwt_mask_fname("c", parameters);
    sdsl::store_to_file(c_mask, fpath);

    fpath = bwt_mask_fname("g", parameters);
    sdsl::store_to_file(g_mask, fpath);

    fpath = bwt_mask_fname("t", parameters);
    sdsl::store_to_file(t_mask, fpath);
}


sdsl::bit_vector load_base_bwt_mask(const std::string &base_char,
                                    const Parameters &parameters) {
    auto fpath = bwt_mask_fname(base_char, parameters);
    sdsl::bit_vector mask;
    sdsl::load_from_file(mask, fpath);
    return mask;
}


DNA_BWT_Masks load_dna_bwt_masks(const FM_Index &fm_index,
                                 const Parameters &parameters) {
    DNA_BWT_Masks dna_bwt_masks;
    dna_bwt_masks.mask_a = load_base_bwt_mask("a", parameters);
    dna_bwt_masks.mask_c = load_base_bwt_mask("c", parameters);
    dna_bwt_masks.mask_g = load_base_bwt_mask("g", parameters);
    dna_bwt_masks.mask_t = load_base_bwt_mask("t", parameters);
    return dna_bwt_masks;
}
