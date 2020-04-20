#include <prg/linearised_prg.hpp>
#include <common/utils.hpp>


/**********************
 * Supporting nesting**
 **********************/
PRG_String::PRG_String(std::string const &file_in, endianness en) :
    odd_site_end_found(false), en(en) {
    std::fstream input(file_in, std::ios::in | std::ios::binary);
    if (!input) throw std::ios::failure("PRG String file not found");
    else {
        uint32_t c{0};
        uint8_t buffer[gram::num_bytes_per_integer]; // Note, use uint8_t and not char which is signed
        size_t byte_pos;
        while (true) {
            input.read((char *) &buffer, gram::num_bytes_per_integer); // Read byte by byte
            if (input.eof()) break;

            if (en == endianness::big) {
                c = (uint32_t) buffer[0] << 24 | (uint32_t) buffer[1] << 16
                    | (uint32_t) buffer[2] << 8 | (uint32_t) buffer[3];
            } else {
                c = (uint32_t) buffer[0] | (uint32_t) buffer[1] << 8
                    | (uint32_t) buffer[2] << 16 | (uint32_t) buffer[3] << 24;
            }

            if (input.eof()) break;
            assert(c >= 1);
            my_PRG_string.push_back(c);
            c = 0;
        }
    };
    output_file = file_in;
    map_and_normalise_ends();

    // Rewrite `file_in`, in two cases:
    // - The PRG string was in legacy format (5G6C5)
    // - sdsl requires it in little endian for building the fm_index
    if (odd_site_end_found || en == endianness::big) write(output_file, endianness::little);
}

PRG_String::PRG_String(marker_vec const &v_in) {
    my_PRG_string = std::move(v_in);
    map_and_normalise_ends();
};

void PRG_String::map_and_normalise_ends() {
    int pos = 0;
    std::size_t v_size = my_PRG_string.size(); // Converts to signed
    Marker marker;
    std::set<Marker> seen_sites;

    while (pos < v_size) {
        marker = my_PRG_string[pos];
        if (marker <= 4) {
            ++pos;
            continue;
        }
        if (is_site_marker(marker)) {
            bool seen_before = seen_sites.find(marker) != seen_sites.end();
            if (seen_before) {
                odd_site_end_found = true;
                end_positions[marker + 1] = pos;
                my_PRG_string[pos]++; // Convert the odd boundary to an even one.
            } else seen_sites.insert(marker);
        } else {
            end_positions[marker] = pos; // Inserts if does not exist, updates otherwise
        }
        ++pos;
    }
}

void PRG_String::write(std::string const &fname, endianness en) {
    std::ofstream out(fname, std::ofstream::out | std::ofstream::binary);
    if (!out) {
        std::cerr << "Cannot open file: " << fname << std::endl;
        exit(1);
    }

    std::vector<char> buffer(gram::num_bytes_per_integer);
    uint8_t byte_number;

    uint8_t model_byte_number;
    if (en == endianness::big) model_byte_number = gram::num_bytes_per_integer - 1;
    else model_byte_number = 0;

    for (auto &s : my_PRG_string) {
        buffer.clear();
        byte_number = model_byte_number;
        while (true) {
            // In big endian, will push the most significant bytes first in buffer
            // In little ^^,  will push the least ^^         ^^      ^^
            buffer.push_back(s >> 8 * byte_number);
            if (en == endianness::big) {
                if (byte_number-- == 0) break;
            } else {
                if (byte_number++ == gram::num_bytes_per_integer - 1) break;
            }
        }
        for (auto e : buffer) out.write(&e, 1);
    }
    out.close();
}

bool operator==(PRG_String const &first, PRG_String const &second) {
    auto const &p_1 = first.get_PRG_string();
    auto const &p_2 = second.get_PRG_string();
    if (p_1.size() != p_2.size()) return false;

    for (int i = 0; i < p_1.size(); ++i) {
        if (p_1[i] != p_2[i]) return false;
    }

    return true;
}

/**************************
 * Not Supporting nesting**
 **************************/

/**
 * Converts a sequence of digits (0-9) into a single integer.
 */
uint64_t concat_marker_digits(const std::vector<int> &marker_digits) {
    uint64_t marker = 0;
    for (const auto &digit: marker_digits)
        marker = marker * 10 + digit;
    return marker;
}

/**
 * Write out marker digits to the encoded prg as a single integer.
 * @see concat_marker_digits()
 */
void flush_marker_digits(std::vector<int> &marker_digits,
                         marker_vec &encoded_prg,
                         uint64_t &count_chars) {
    if (marker_digits.empty())
        return;

    uint64_t marker = concat_marker_digits(marker_digits);
    encoded_prg[count_chars++] = marker;
    marker_digits.clear();
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

std::string gram::ints_to_prg_string(std::vector<Marker> const& int_vec){
    std::string readable_string(int_vec.size(), '0');
    std::unordered_map<int,int> last_allele_indices; // Will record where to close the sites.

    int pos{-1};
    for (auto& s : int_vec){
        pos++;
        if (s > 4){
            if (s%2 == 1) readable_string[pos] = '[';
            else{
                readable_string[pos] = ',';
                if (last_allele_indices.find(s) != last_allele_indices.end()){
                    last_allele_indices.erase(s);
                }
                last_allele_indices.insert({s, pos});
            }
            continue;
        }
        // Implicit else: s <= 4
        char base = decode_dna_base(s).c_str()[0];
        readable_string[pos] = base;
    }

    // Close the sites
    for (auto s : last_allele_indices){
        auto pos = s.second;
        readable_string[pos] = ']';
    }
    return readable_string;
}

marker_vec gram::prg_string_to_ints(std::string const& string_prg) {
    int_Base base;
    std::stack<int> marker_stack;
    int max_var_marker{3};
    int char_count{0};

    marker_vec encoded_prg(string_prg.size());
    for (std::size_t i = 0; i < string_prg.size(); ++i){
        const auto &c = string_prg[i];

        switch(c) {
            case '[' : {
                max_var_marker += 2;
                marker_stack.push(max_var_marker);
                encoded_prg[char_count++] = max_var_marker;
                break;
            }

            case ']' : {
                assert(!marker_stack.empty());
                encoded_prg[char_count++] = marker_stack.top() + 1;
                marker_stack.pop();
                break;
            }

            case ',' : {
                assert(!marker_stack.empty());
                encoded_prg[char_count++] = marker_stack.top() + 1;
                break;
            }

            default : {
                base = encode_dna_base(c);
                if (base == 0){
                    std::cerr << "Error: the argument " << c << " is not a nucleotide char";
                    exit(1);
                }
                encoded_prg[char_count++] = encode_dna_base(c);
                break;
            }
        }
    }

    // BOOST_LOG_TRIVIAL(info) << "Number of sites produced: " << (max_var_marker -3 ) / 2;
    return encoded_prg;
}
