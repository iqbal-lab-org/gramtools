#include <prg/linearised_prg.hpp>


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
};

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
