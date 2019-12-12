#include "prg/masks.hpp"
#include "prg/prg.hpp"


using namespace gram;

uint64_t gram::dna_bwt_rank(const uint64_t &upper_index,
                            const Marker &dna_base,
                            const PRG_Info &prg_info) {
    switch (dna_base) {
        case 1:
            return prg_info.rank_bwt_a(upper_index);
        case 2:
            return prg_info.rank_bwt_c(upper_index);
        case 3:
            return prg_info.rank_bwt_g(upper_index);
        case 4:
            return prg_info.rank_bwt_t(upper_index);
        default:
            return 0;
    }
}


//marker_vec gram::generate_encoded_prg(const Parameters &parameters) {
//    auto encoded_prg = parse_raw_prg_file(parameters.linear_prg_fpath);
//    sdsl::store_to_file(encoded_prg, parameters.encoded_prg_fpath);
//    return encoded_prg;
//}


marker_vec gram::parse_raw_prg_file(const std::string &prg_fpath) {
    const auto prg_raw = load_raw_prg(prg_fpath);
    auto encoded_prg = encode_prg(prg_raw);
    return encoded_prg;
}

std::string gram::load_raw_prg(const std::string &prg_fpath) {
    std::ifstream fhandle(prg_fpath, std::ios::in | std::ios::binary);
    if (not fhandle) {
        std::cout << "Problem reading PRG input file" << std::endl;
        exit(1);
    }

    std::string prg;
    // Sets the input position indicator to end of file (0 relative to end).
    fhandle.seekg(0, std::ios::end);
    // Set the prg string size to the number of characters in the file.
    prg.resize((unsigned long) fhandle.tellg());

    // Go back to file beginning, and read all characters to prg string.
    fhandle.seekg(0, std::ios::beg);
    fhandle.read(&prg[0], prg.size());
    fhandle.close();

    return prg;
}

marker_vec gram::encode_prg(const std::string &prg_raw) {
    marker_vec encoded_prg(prg_raw.length(), 0);

    uint64_t count_chars = 0;
    // TODO: this should be possible without storing each individual digit
    std::vector<int> marker_digits;
    for (const auto &c: prg_raw) {
        EncodeResult encode_result = encode_char(c);

        if (encode_result.is_dna) {
            // `flush_marker_digits` flushes any latent marker characters
            flush_marker_digits(marker_digits, encoded_prg, count_chars);
            encoded_prg[count_chars++] = encode_result.character;
            continue;
        }

        // else: record the digit, and stand ready to record another
        // TODO: check that character is numeric?
        marker_digits.push_back(encode_result.character);
    }
    flush_marker_digits(marker_digits, encoded_prg, count_chars);

    encoded_prg.resize(count_chars);
    return encoded_prg;
}

void gram::flush_marker_digits(std::vector<int> &marker_digits,
                               marker_vec &encoded_prg,
                               uint64_t &count_chars) {
    if (marker_digits.empty())
        return;

    uint64_t marker = concat_marker_digits(marker_digits);
    encoded_prg[count_chars++] = marker;
    marker_digits.clear();
}

uint64_t gram::concat_marker_digits(const std::vector<int> &marker_digits) {
    uint64_t marker = 0;
    for (const auto &digit: marker_digits)
        marker = marker * 10 + digit;
    return marker;
}

EncodeResult gram::encode_char(const char &c) {
    EncodeResult encode_result = {};

    switch (c) {
        case 'A':
        case 'a':
            encode_result.is_dna = true;
            encode_result.character = 1;
            return encode_result;

        case 'C':
        case 'c':
            encode_result.is_dna = true;
            encode_result.character = 2;
            return encode_result;

        case 'G':
        case 'g':
            encode_result.is_dna = true;
            encode_result.character = 3;
            return encode_result;

        case 'T':
        case 't':
            encode_result.is_dna = true;
            encode_result.character = 4;
            return encode_result;

        default:
            /// The character is non-DNA, so must be a variant marker.
            encode_result.is_dna = false;
            encode_result.character = (uint32_t) c - '0';
            return encode_result;
    }
}

PRG_Info gram::load_prg_info(const Parameters &parameters) {
    PRG_Info prg_info;

    PRG_String ps{parameters.encoded_prg_fpath};
    prg_info.last_allele_positions = ps.get_end_positions();

    // Load coverage graph
    std::ifstream ifs{parameters.cov_graph_fpath};
    boost::archive::binary_iarchive ia{ifs};
    ia >> prg_info.coverage_graph;
    prg_info.num_variant_sites = prg_info.coverage_graph.bubble_map.size();

    prg_info.fm_index = load_fm_index(parameters);

    prg_info.bwt_markers_mask = generate_bwt_markers_mask(prg_info.fm_index);

    prg_info.dna_bwt_masks = load_dna_bwt_masks(prg_info.fm_index, parameters);
    prg_info.rank_bwt_a = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_a);
    prg_info.rank_bwt_c = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_c);
    prg_info.rank_bwt_g = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_g);
    prg_info.rank_bwt_t = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_t);

    return prg_info;
}