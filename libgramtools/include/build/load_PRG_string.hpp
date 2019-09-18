#include <string>
#include <set>
#include <fstream>
#include <cstring>
#include <iostream>
#include <assert.h>
#include "common/utils.hpp"


using namespace gram;


class PRG_String {
public:
    PRG_String() {};

    /**
     * Read in PRG String from binary int vector
     * Reads byte per byte in BIG ENDIAN; the serialisor must write that way too.
     */
    PRG_String(std::string const &file_in);

    PRG_String(marker_vec const &v_in);

    /*
     * Functions
     */

    /**
     * PRG string traversal:
     * -first to figure out the site ends
     * -second to get parental relationships
     */
    void process();

    void write();

    void set_output_file(std::string const &fname) { output_file = fname; };

    // Getters
    const marker_vec get_PRG_string() const { return my_PRG_string; };
    const std::unordered_map<Marker, int> get_end_positions() const { return end_positions;};

    friend std::ostream &operator<<(std::ostream &out, PRG_String const &e) {
        for (auto &s : e.my_PRG_string) out << s;
        return out;
    }

    /*
     * Variables
     */
    bool odd_site_end_found; // If the case, we will rewrite the int vector.

private:
    std::string output_file;
    marker_vec my_PRG_string;
    std::unordered_map<Marker, int> end_positions; // Where a given site ends; the allele (even) marker is stored.

    /**
     * Discover where site boundaries lie, and convert any odd end markers to even end markers
     */
    void map_and_normalise_ends();
};
