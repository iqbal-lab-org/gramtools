#include <boost/filesystem.hpp>

#include "common/utils.hpp"
#include "prg/dna_ranks.hpp"


namespace fs = boost::filesystem;
using namespace gram;

/**
 * Generates a bit vector with bit set if the given DNA `base` is present at each index of the BWT.
 */
void populate_dna_bwt_masks(FM_Index const& fm_index, DNA_BWT_Masks& d_m) {
    for (uint64_t i = 0; i < fm_index.bwt.size(); i++){
        switch(fm_index.bwt[i]){
            case 1:
                d_m.mask_a[i] = true;
                break;
            case 2:
                d_m.mask_c[i] = true;
                break;
            case 3:
                d_m.mask_g[i] = true;
                break;
            case 4:
                d_m.mask_t[i] = true;
                break;
        }
    }
}

/**
 * Generates a filename for a BWT mask for nucleotide bases.
 * @see generate_base_bwt_mask()
 */
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

DNA_BWT_Masks gram::generate_bwt_masks(FM_Index const& fm_index,
                              const Parameters &parameters) {
    auto bwt_size = fm_index.bwt.size();
    DNA_BWT_Masks d_m;
    // Set storage for the bit masks
    d_m.mask_a = sdsl::bit_vector(bwt_size, 0);
    d_m.mask_c = sdsl::bit_vector(bwt_size, 0);
    d_m.mask_g = sdsl::bit_vector(bwt_size, 0);
    d_m.mask_t = sdsl::bit_vector(bwt_size, 0);

    populate_dna_bwt_masks(fm_index, d_m);


    auto fpath = bwt_mask_fname("a", parameters);
    sdsl::store_to_file(d_m.mask_a, fpath);

    fpath = bwt_mask_fname("c", parameters);
    sdsl::store_to_file(d_m.mask_c, fpath);

    fpath = bwt_mask_fname("g", parameters);
    sdsl::store_to_file(d_m.mask_g, fpath);

    fpath = bwt_mask_fname("t", parameters);
    sdsl::store_to_file(d_m.mask_t, fpath);

    return d_m;
}


sdsl::bit_vector load_base_bwt_mask(const std::string &base_char,
                                    const Parameters &parameters) {
    auto fpath = bwt_mask_fname(base_char, parameters);
    sdsl::bit_vector mask;
    sdsl::load_from_file(mask, fpath);
    return mask;
}


DNA_BWT_Masks gram::load_dna_bwt_masks(const FM_Index &fm_index,
                                       const Parameters &parameters) {
    DNA_BWT_Masks dna_bwt_masks;
    dna_bwt_masks.mask_a = load_base_bwt_mask("a", parameters);
    dna_bwt_masks.mask_c = load_base_bwt_mask("c", parameters);
    dna_bwt_masks.mask_g = load_base_bwt_mask("g", parameters);
    dna_bwt_masks.mask_t = load_base_bwt_mask("t", parameters);
    return dna_bwt_masks;
}
