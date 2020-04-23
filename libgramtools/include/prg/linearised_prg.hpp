/**
 * @file
 * In gramtools, we use a linearised representation of a Population Reference Graph (PRG) for mapping reads to.
 * The one that supports arbitrarily nested variation reads such linear prgs from a stream of binary integers.
 */
#ifndef PRG_STRING_HPP
#define PRG_STRING_HPP

#include <string>
#include <set>
#include <fstream>
#include <cstring>
#include <iostream>
#include <assert.h>

#include "common/data_types.hpp"

using namespace gram;

namespace gram{
    enum class endianness{big, little};

    /** Converts linearised PRG as int vector to a more readable string.
     *  We use the following notation: '[' opens a site, ',' delimits alleles in a site, ']' closes a site.
     * */
    std::string ints_to_prg_string(std::vector<Marker> const& int_vec);


    /**
     * Convert a nested PRG string to int representation, with linear site numbering.
     * The site numbering is based on the fixed order in which '[' chars are encountered;
     * thus int -> prg_string -> int can lose original site numbering.
     */
    std::vector<Marker> prg_string_to_ints(std::string const& string_prg);
}

/**********************
 * Supporting nesting**
 **********************/
class PRG_String {
public:
    PRG_String() = default;

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
    marker_vec get_PRG_string() const { return my_PRG_string; };
    std::size_t size() const { return my_PRG_string.size(); };
    endianness get_endianness() const { return en; };
    std::unordered_map<Marker, int> get_end_positions() const { return end_positions;};

    friend std::ostream &operator<<(std::ostream &out, PRG_String const &e) {
        for (auto &s : e.my_PRG_string) out << s;
        return out;
    }

    friend bool operator==(PRG_String const& first, PRG_String const& second);

    /*
     * Variables
     */
    bool odd_site_end_found{false}; // If the case, we will rewrite the int vector.

private:
    std::string output_file;
    endianness en{endianness::little};
    marker_vec my_PRG_string;
    std::unordered_map<Marker, int> end_positions; // Where a given site ends; the allele (even) marker is stored.

    /**
     * Discover where site boundaries lie, and convert any odd end markers to even end markers
     */
    void map_ends_and_check_for_duplicates();
};


/**************************
 * Not Supporting nesting**
 **************************/
namespace gram{
/**
 * Convert prg as string of characters to vector of integers.
 * Nucleotides encoded as 1-4. Variant markers can make up several characters so are treated with a buffer.
 * NB: this function only works for PRGs with no nested variation (otherwise, for eg, '57' confounded with '5' then '7')
 */
    marker_vec encode_prg(const std::string &prg_raw);
}

#endif //PRG_STRING_HPP
