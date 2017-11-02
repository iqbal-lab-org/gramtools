#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "parameters.hpp"


#ifndef GRAMTOOLS_PROCESS_PRG_HPP
#define GRAMTOOLS_PROCESS_PRG_HPP

using EncodedPRG = sdsl::int_vector<>;
using WavletTree = sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>>;
using FM_Index = sdsl::csa_wt<WavletTree, 2, 16777216>;

EncodedPRG generate_encoded_prg(const Parameters &parameters);

FM_Index load_fm_index(const Parameters &parameters);

void dump_encoded_prg(const std::vector<uint64_t> &prg,
                      const std::string &prg_encoded_fpath);

FM_Index generate_fm_index(const Parameters &parameters);

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

#endif //GRAMTOOLS_PROCESS_PRG_HPP
