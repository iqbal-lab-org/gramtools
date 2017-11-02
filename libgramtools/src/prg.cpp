#include "masks.hpp"
#include "prg.hpp"


uint64_t get_max_alphabet_num(const sdsl::int_vector<> &encoded_prg) {
    uint64_t max_alphabet_num = 0;
    for (const auto &x: encoded_prg) {
        if (x > max_alphabet_num)
            max_alphabet_num = x;
    }
    return max_alphabet_num;
}


sdsl::int_vector<> generate_encoded_prg(const Parameters &parameters) {
    auto encoded_prg = parse_raw_prg_file(parameters.linear_prg_fpath);
    sdsl::store_to_file(encoded_prg, parameters.encoded_prg_fpath);
    return encoded_prg;
}


sdsl::int_vector<> parse_raw_prg_file(const std::string &prg_fpath) {
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


sdsl::int_vector<> encode_prg(const std::string &prg_raw) {
    sdsl::int_vector<> encoded_prg(prg_raw.length(), 0, 64);

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
                         sdsl::int_vector<> &encoded_prg,
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


PRG_Info load_prg_info(const Parameters &parameters) {
    MasksParser masks(parameters.site_mask_fpath,
                      parameters.allele_mask_fpath);
    auto fm_index = load_fm_index(parameters);
    auto encoded_prg = parse_raw_prg_file(parameters.linear_prg_fpath);
    auto markers_mask = generate_markers_mask(encoded_prg);
    auto max_alphabet_num = get_max_alphabet_num(encoded_prg);
    return PRG_Info {
            fm_index,
            encoded_prg,
            masks.sites,
            masks.allele,
            markers_mask,
            max_alphabet_num
    };
}