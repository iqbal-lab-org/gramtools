#include <prg/load_PRG_string.hpp>


PRG_String::PRG_String(std::string const& file_in){
    std::fstream input(file_in, std::ios::in | std::ios::binary);
    if (!input) throw std::ios::failure("PRG String file not found");
    else {
        int32_t c;
        char buffer[gram::num_bytes_per_integer];
        while (true){
            input.read(buffer, 4); // Read 4 bytes in one go
            if (input.eof()) break;
            // Place the bytes in the int32, assuming BIG ENDIAN
            c = buffer[0] << 24 | buffer[1] << 16 | buffer[2] << 8 | buffer[3];
            assert(c >= 1);
            my_PRG_string.push_back(c);
        }
    };
    output_file = file_in;
    map_and_normalise_ends();
}

PRG_String::PRG_String(marker_vec const& v_in){
   my_PRG_string = std::move(v_in);
   map_and_normalise_ends();
};

void PRG_String::map_and_normalise_ends(){
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

void PRG_String::write(std::string const& fname, endianness en){
 std::ofstream out (fname, std::ofstream::out | std::ofstream::binary);
 if(!out){
     std::cerr << "Cannot open file: " << fname << std::endl;
     exit(1);
 }

 std::vector<char> buffer(gram::num_bytes_per_integer);
 uint8_t byte_number;
uint8_t model_byte_number;
en == endianness::big ? model_byte_number = gram::num_bytes_per_integer - 1 : 0;
 for (auto& s : my_PRG_string){
     buffer.clear();
     byte_number = model_byte_number;
     while (true){
            // In big endian, will push the most significant bytes first in buffer
            // In little ^^,  will push the least ^^         ^^      ^^
            buffer.push_back(s >> 8 * byte_number);
            if (en == endianness::big){
                if (byte_number-- == 0) break;
            }
            else {
                if (byte_number++ == gram::num_bytes_per_integer - 1) break;
            }
     }
     for (auto e : buffer) out.write(&e, 1);
 }
 out.close();
}

bool operator==(PRG_String const& first, PRG_String const& second){
    auto const& p_1 = first.get_PRG_string();
    auto const& p_2 = second.get_PRG_string();
    if (p_1.size() != p_2.size()) return false;

    for (int i = 0; i < p_1.size(); ++i){
        if (p_1[i] != p_2[i]) return false;
    }

    return true;
}
