#include <prg/load_PRG_string.hpp>


PRG_String::PRG_String(std::string const& file_in){
    std::fstream input(file_in, std::ios::in | std::ios::binary);
    if (!input) throw std::ios::failure("PRG String file not found");
    else {
        int32_t c;
        char buffer[4];
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
    int pos = -1;
    int v_size = my_PRG_string.size(); // Converts to signed
    Marker marker;
    std::set<Marker> seen_sites;

    while (pos < v_size) {
        marker = my_PRG_string[++pos];
        if (marker <= 4) continue;
        if (marker % 2 == 1) {
            bool seen_before = seen_sites.find(marker) != seen_sites.end();
            if (seen_before) {
                odd_site_end_found = true;
                end_positions.insert(std::make_pair(marker + 1, pos));
                my_PRG_string[pos]++; // Convert the odd boundary to an even one.
            } else seen_sites.insert(marker);
        } else {
            try {
                end_positions.at(marker) = pos;
            }
            catch (std::out_of_range &c) {
                end_positions.insert(std::make_pair(marker, pos));
            }
        }
    }
};
