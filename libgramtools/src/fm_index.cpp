#include <fstream>
#include <cinttypes>
#include <sdsl/suffix_arrays.hpp>

#include "parameters.hpp"
#include "fm_index.hpp"


EncodedPRG generate_encoded_prg(const Parameters &parameters) {
    auto encoded_prg = parse_prg(parameters.linear_prg_fpath);
    sdsl::store_to_file(encoded_prg, parameters.encoded_prg_fpath);
    return encoded_prg;
}


FM_Index load_fm_index(const Parameters &parameters) {
    FM_Index fm_index;
    sdsl::load_from_file(fm_index, parameters.fm_index_fpath);
    return fm_index;
}


FM_Index generate_fm_index(const Parameters &parameters) {
    FM_Index fm_index;

    sdsl::memory_monitor::start();
    sdsl::construct(fm_index, parameters.encoded_prg_fpath, 0);
    sdsl::memory_monitor::stop();

    std::ofstream memory_log_fhandle(parameters.sdsl_memory_log_fpath);
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(memory_log_fhandle);

    sdsl::store_to_file(fm_index, parameters.fm_index_fpath);
    return fm_index;
}


EncodedPRG parse_prg(const std::string &prg_fpath) {
    const auto prg_raw = load_raw_prg(prg_fpath);
    auto encoded_prg = encode_prg(prg_raw);
    return encoded_prg;
}


std::string load_raw_prg(const std::string &prg_fpath) {
    std::ifstream fhandle(prg_fpath, std::ios::in | std::ios::binary);
    if (not fhandle) {
        std::cout << "Problem reading PRG input file" << std::endl;
        exit(1);
    }

    std::string prg;
    fhandle.seekg(0, std::ios::end);
    prg.resize((unsigned long) fhandle.tellg());

    fhandle.seekg(0, std::ios::beg);
    fhandle.read(&prg[0], prg.size());
    fhandle.close();

    return prg;
}


EncodedPRG encode_prg(const std::string &prg_raw) {
    EncodedPRG encoded_prg(prg_raw.length(), 0, 64);

    uint64_t count_chars = 0;
    // TODO: this should be possible without storing each individual digit
    std::vector<int> marker_digits;
    for (const auto &c: prg_raw) {
        EncodeResult encode_result = encode_char(c);

        if (encode_result.is_dna) {
            flush_marker_digits(marker_digits, encoded_prg, count_chars);
            encoded_prg[count_chars++] = encode_result.charecter;
            continue;
        }

        marker_digits.push_back(encode_result.charecter);
    }

    flush_marker_digits(marker_digits, encoded_prg, count_chars);
    encoded_prg.resize(count_chars);
    return encoded_prg;
}


void flush_marker_digits(std::vector<int> &marker_digits,
                         EncodedPRG &encoded_prg,
                         uint64_t &count_chars) {
    if (marker_digits.empty())
        return;

    uint64_t marker = concat_marker_digits(marker_digits);
    encoded_prg[count_chars++] = marker;
    marker_digits.clear();
}


uint64_t concat_marker_digits(const std::vector<int> &marker_digits) {
    uint64_t marker = 0;
    for (const auto &digit: marker_digits)
        marker = marker * 10 + digit;
    return marker;
}


EncodeResult encode_char(const char &c) {
    EncodeResult encode_result;

    switch (c) {
        case 'A':
        case 'a':
            encode_result.is_dna = true;
            encode_result.charecter = 1;
            return encode_result;

        case 'C':
        case 'c':
            encode_result.is_dna = true;
            encode_result.charecter = 2;
            return encode_result;

        case 'G':
        case 'g':
            encode_result.is_dna = true;
            encode_result.charecter = 3;
            return encode_result;

        case 'T':
        case 't':
            encode_result.is_dna = true;
            encode_result.charecter = 4;
            return encode_result;

        default:
            encode_result.is_dna = false;
            encode_result.charecter = c - '0';
            return encode_result;
    }
}
