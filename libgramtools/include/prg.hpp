#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>
#include <tuple>

#include "utils.hpp"
#include "fm_index.hpp"


#ifndef GRAMTOOLS_PRG_HPP
#define GRAMTOOLS_PRG_HPP

struct PRG_Info {
    FM_Index fm_index;
    sdsl::int_vector<> encoded_prg;
    // DNA_Rank dna_rank;
    std::vector<Marker> sites_mask;
    std::vector<AlleleId> allele_mask;
    sdsl::bit_vector markers_mask;
    uint64_t max_alphabet_num;
};

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

PRG_Info load_prg_info(const Parameters &params);

#endif //GRAMTOOLS_PRG_HPP
