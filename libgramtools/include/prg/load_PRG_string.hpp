#ifndef PRG_STRING_HPP
#define PRG_STRING_HPP
#include <string>
#include <set>
#include <fstream>
#include <cstring>
#include <iostream>
#include <assert.h>
#include "common/utils.hpp"
#include "common/parameters.hpp"


using namespace gram;

namespace gram{
    enum class endianness{big, little};
}

class PRG_String {
public:
    PRG_String() {};

    /**
     * Read in PRG String from binary int vector
     * Reads byte per byte in specified endianness; the serialisor must write that way too.
     */
    PRG_String(std::string const &file_in, endianness en = endianness::little);

    PRG_String(marker_vec const &v_in);

    /*
     * Functions
     */
    void write(std::string const& fname, endianness = endianness::little);

    // Getters
    const marker_vec get_PRG_string() const { return my_PRG_string; };
    const endianness get_endianness() const { return en; };
    const std::unordered_map<Marker, int> get_end_positions() const { return end_positions;};

    friend std::ostream &operator<<(std::ostream &out, PRG_String const &e) {
        for (auto &s : e.my_PRG_string) out << s;
        return out;
    }

    friend bool operator==(PRG_String const& first, PRG_String const& second);

    /*
     * Variables
     */
    bool odd_site_end_found; // If the case, we will rewrite the int vector.

private:
    std::string output_file;
    endianness en;
    marker_vec my_PRG_string;
    std::unordered_map<Marker, int> end_positions; // Where a given site ends; the allele (even) marker is stored.

    /**
     * Discover where site boundaries lie, and convert any odd end markers to even end markers
     */
    void map_and_normalise_ends();
};
#endif //COV_GRAPH_HPP
