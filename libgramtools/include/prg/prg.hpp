#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>
#include <tuple>

#include "common/utils.hpp"
#include "dna_ranks.hpp"
#include "fm_index.hpp"


#ifndef GRAMTOOLS_PRG_HPP
#define GRAMTOOLS_PRG_HPP


struct PRG_Info {
    FM_Index fm_index;
    sdsl::int_vector<> encoded_prg;
    
    std::vector<Marker> sites_mask;
    sdsl::int_vector<> allele_mask;

    sdsl::bit_vector bwt_markers_mask;
    sdsl::rank_support_v<1> bwt_markers_rank;
    sdsl::select_support_mcl<1> bwt_markers_select;
    uint64_t markers_mask_count_set_bits;

    sdsl::bit_vector prg_markers_mask;
    sdsl::rank_support_v<1> prg_markers_rank;
    sdsl::select_support_mcl<1> prg_markers_select;

    DNA_BWT_Masks dna_bwt_masks;
    sdsl::rank_support_v<1> rank_bwt_a;
    sdsl::rank_support_v<1> rank_bwt_c;
    sdsl::rank_support_v<1> rank_bwt_g;
    sdsl::rank_support_v<1> rank_bwt_t;

    uint64_t max_alphabet_num;
};

uint64_t dna_bwt_rank(const uint64_t upper_index,
                      const Marker &dna_base,
                      const PRG_Info &prg_info);

uint64_t get_max_alphabet_num(const sdsl::int_vector<> &encoded_prg);

sdsl::int_vector<> generate_encoded_prg(const Parameters &parameters);

sdsl::int_vector<> parse_raw_prg_file(const std::string &prg_fpath);

std::string load_raw_prg(const std::string &prg_fpath);

sdsl::int_vector<> encode_prg(const std::string &prg_raw);

void flush_marker_digits(std::vector<int> &marker_digits,
                         sdsl::int_vector<> &encoded_prg,
                         uint64_t &count_chars);

uint64_t concat_marker_digits(const std::vector<int> &marker_digits);

struct EncodeResult{
    bool is_dna;
    int charecter;
};

EncodeResult encode_char(const char &c);

PRG_Info load_prg_info(const Parameters &parameters);

#endif //GRAMTOOLS_PRG_HPP
