#include <fstream>
#include <cinttypes>
#include <sdsl/suffix_arrays.hpp>

#include "parameters.hpp"
#include "fm_index.hpp"


void generate_encoded_prg(const Parameters &parameters) {
    auto encoded_prg = parse_prg(parameters.linear_prg_fpath);
    std::cout << "Number of charecters in integer encoded linear PRG: "
              << encoded_prg.size() << std::endl;
    dump_encoded_prg(encoded_prg, parameters.encoded_prg_fpath);
}


void dump_encoded_prg(const EncodedPRG &prg,
                      const std::string &prg_encoded_fpath) {
    std::ofstream fhandle;
    fhandle.open(prg_encoded_fpath, std::ios::out | std::ios::binary);
    fhandle.write((char *) &prg[0], prg.size() * sizeof(PRG_Char));
    fhandle.close();
}


FM_Index load_fm_index(const Parameters &parameters) {
    FM_Index fm_index;
    sdsl::load_from_file(fm_index, parameters.fm_index_fpath);
    return fm_index;
}


void generate_fm_index(const Parameters &parameters) {
    FM_Index fm_index;

    sdsl::memory_monitor::start();
    sdsl::construct(fm_index, parameters.encoded_prg_fpath, 8);
    sdsl::memory_monitor::stop();

    std::ofstream memory_log_fhandle(parameters.sdsl_memory_log_fpath);
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(memory_log_fhandle);

    sdsl::store_to_file(fm_index, parameters.fm_index_fpath);
}


EncodedPRG parse_prg(const std::string &prg_fpath) {
    const auto prg_raw = load_raw_prg(prg_fpath);
    auto prg = encode_prg(prg_raw);
    return prg;
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
    EncodedPRG prg_encoded;
    prg_encoded.reserve(prg_raw.length());

    // TODO: this should be possible without storing each individual digit
    std::vector<int> marker_digits;
    for (const auto &c: prg_raw) {
        EncodeResult encode_result = encode_char(c);

        if (encode_result.is_dna) {
            flush_marker_digits(marker_digits, prg_encoded);
            prg_encoded.emplace_back(encode_result.charecter);
            continue;
        }

        marker_digits.push_back(encode_result.charecter);
    }

    flush_marker_digits(marker_digits, prg_encoded);
    prg_encoded.shrink_to_fit();
    return prg_encoded;
}


void flush_marker_digits(std::vector<int> &marker_digits,
                         EncodedPRG &encoded_prg) {
    if (marker_digits.empty())
        return;

    uint64_t marker = concat_marker_digits(marker_digits);
    encoded_prg.push_back(marker);
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
