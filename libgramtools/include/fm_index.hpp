//
// Created by robyn on 09/06/17.
//

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#ifndef GRAMTOOLS_PROCESS_PRG_HPP
#define GRAMTOOLS_PROCESS_PRG_HPP


using WavletTree = sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>>;
using FM_Index = sdsl::csa_wt<WavletTree, 2, 16777216>;

FM_Index construct_fm_index(bool fwd, std::string fm_index_fpath, std::string prg_encoded_fpath, const std::string &prg_fpath,
                            const std::string &memory_log_fname);

void dump_encoded_prg(const std::vector<uint64_t> &prg,
                      const std::string &prg_encoded_fpath);

FM_Index build_fm_index(const std::string &prg_encoded_fpath,
                        const std::string &fm_index_fpath,
                        const std::string &memory_log_fname);

std::vector<uint64_t> parse_prg(const std::string &prg_fpath);

std::string load_raw_prg(const std::string &prg_fpath);

std::vector<uint64_t> encode_prg(const std::string &prg_raw);

void flush_marker_digits(std::vector<int> &marker_digits, std::vector<uint64_t> &prg_encoded);

uint64_t concat_marker_digits(const std::vector<int> &marker_digits);

struct EncodeResult{
    bool is_dna;
    int charecter;
};

EncodeResult encode_char(const char &c);


#endif //GRAMTOOLS_PROCESS_PRG_HPP
