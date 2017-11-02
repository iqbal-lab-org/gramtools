#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>
#include <tuple>

#include "masks.hpp"
#include "fm_index.hpp"
#include "ranks.hpp"


#ifndef GRAMTOOLS_PRG_HPP
#define GRAMTOOLS_PRG_HPP

using EncodedPRG = sdsl::int_vector<>;

struct PRG_Info {
    FM_Index fm_index;
    EncodedPRG encoded_prg;
    // DNA_Rank dna_rank;
    SitesMask sites_mask;
    AlleleMask allele_mask;
    uint64_t max_alphabet_num;
};

uint64_t get_max_alphabet_num(const EncodedPRG &encoded_prg);

EncodedPRG generate_encoded_prg(const Parameters &parameters);

EncodedPRG parse_prg(const std::string &prg_fpath);

std::string load_raw_prg(const std::string &prg_fpath);

EncodedPRG encode_prg(const std::string &prg_raw);

void flush_marker_digits(std::vector<int> &marker_digits,
                         EncodedPRG &encoded_prg,
                         uint64_t &count_chars);

uint64_t concat_marker_digits(const std::vector<int> &marker_digits);

struct EncodeResult{
    bool is_dna;
    int charecter;
};

EncodeResult encode_char(const char &c);

PRG_Info load_prg_info(const Parameters &params);

#endif //GRAMTOOLS_PRG_HPP
