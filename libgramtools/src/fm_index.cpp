#include <fstream>
#include <cinttypes>
#include <sdsl/suffix_arrays.hpp>

#include "fm_index.hpp"
#include "map.hpp"


//make SA sampling density and ISA sampling density customizable
//make void fcn and pass csa by reference? return ii?
FM_Index construct_fm_index(bool fwd,
                            std::string fm_index_fpath,
                            std::string prg_encoded_fpath,
                            const std::string &prg_fpath,
                            const std::string &memory_log_fname) {

    std::vector<uint64_t> prg = parse_prg(prg_fpath);
    std::cout << "Number of integers in int encoded linear PRG: " << prg.size() << std::endl;

    if (!fwd) {
        prg_encoded_fpath = prg_encoded_fpath + "_rev";
        fm_index_fpath = fm_index_fpath + "_rev";
        std::reverse(prg.begin(), prg.end());
    }

    dump_encoded_prg(prg, prg_encoded_fpath);
    FM_Index fm_index = build_fm_index(prg_encoded_fpath,
                                       fm_index_fpath,
                                       memory_log_fname);
    return fm_index;
}


void dump_encoded_prg(const std::vector<uint64_t> &prg,
                      const std::string &prg_encoded_fpath) {
    std::ofstream fhandle;
    fhandle.open(prg_encoded_fpath, std::ios::out | std::ios::binary);
    fhandle.write((char *) &prg[0], prg.size() * sizeof(uint64_t));
    fhandle.close();
}


bool file_exist(const std::string &fpath) {
    std::ifstream fhandle(fpath.c_str());
    auto result = fhandle.good();
    fhandle.close();
    return result;
}


FM_Index build_fm_index(const std::string &prg_encoded_fpath,
                        const std::string &fm_index_fpath,
                        const std::string &memory_log_fname) {

    FM_Index fm_index;

    sdsl::memory_monitor::start();
    sdsl::construct(fm_index, prg_encoded_fpath, 8);
    sdsl::memory_monitor::stop();

    std::ofstream memory_log_fhandle(memory_log_fname);
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(memory_log_fhandle);

    sdsl::store_to_file(fm_index, fm_index_fpath);

    assert(!fm_index.empty());
    return fm_index;

    /*

    if(file_exist(fm_index_fpath)) {
        load_from_file(fm_index, fm_index_fpath);
        if (!fm_index.empty())
            return fm_index;
    }

    sdsl::memory_monitor::start();
    sdsl::construct(fm_index, prg_encoded_fpath, 8);
    sdsl::memory_monitor::stop();

    std::ofstream memory_log_fhandle(memory_log_fname);
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(memory_log_fhandle);

    sdsl::store_to_file(fm_index, fm_index_fpath);

    assert(!fm_index.empty());
    return fm_index;
     */
}


std::vector<uint64_t> parse_prg(const std::string &prg_fpath) {
    std::string prg_raw = load_raw_prg(prg_fpath);
    std::vector<uint64_t> prg = encode_prg(prg_raw);
    return prg;
}


std::string load_raw_prg(const std::string &prg_fpath) {
    std::ifstream fhandle(prg_fpath, std::ios::in | std::ios::binary);
    if (!fhandle) {
        std::cout << "Problem reading PRG input file" << std::endl;
        exit(1);
    }

    std::string prg;
    fhandle.seekg(0, std::ios::end);
    prg.resize(fhandle.tellg());

    fhandle.seekg(0, std::ios::beg);
    fhandle.read(&prg[0], prg.size());
    fhandle.close();

    return prg;
}


std::vector<uint64_t> encode_prg(const std::string &prg_raw) {
    std::vector<uint64_t> prg_encoded;
    prg_encoded.reserve(prg_raw.length());

    // TODO: this should be possible without storing each individual digit
    std::vector<int> marker_digits;
    for (const auto &c: prg_raw) {
        EncodeResult encode_result = encode_char(c);

        if (encode_result.is_dna) {
            flush_marker_digits(marker_digits, prg_encoded);
            prg_encoded.push_back(encode_result.charecter);
            continue;
        }

        marker_digits.push_back(encode_result.charecter);
    }

    flush_marker_digits(marker_digits, prg_encoded);
    prg_encoded.shrink_to_fit();
    return prg_encoded;
}


void flush_marker_digits(std::vector<int> &marker_digits,
                         std::vector<uint64_t> &prg_encoded) {
    if (marker_digits.empty())
        return;

    uint64_t marker = concat_marker_digits(marker_digits);
    prg_encoded.push_back(marker);
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
