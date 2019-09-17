#include <unordered_map>
#include <string>
#include <set>
#include <vector>
#include <fstream>
#include <cstring>
#include <iostream>
#include <assert.h>

using par_map = std::unordered_map<int, std::pair<int,int>>;
using Marker = int32_t;
using int_vec = std::vector<Marker>;

class PRG_String{
public:
    PRG_String(){};
    /**
     * Read in PRG String from binary int vector
     * Reads byte per byte in BIG ENDIAN; the serialisor must write that way too.
     */
    PRG_String(std::string const& file_in);
    PRG_String(int_vec const& v_in);

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
    void set_output_file(std::string const& fname) {output_file = fname;};
    const int_vec get_PRG_string() const {return my_PRG_string;};

    friend std::ostream& operator<<(std::ostream& out, PRG_String const& e){
        for (auto&s : e.my_PRG_string) out << s;
        return out;
    }

    /*
     * Variables
     */
     bool odd_site_end_found; // If the case, we will rewrite the int vector.
    std::unordered_map<Marker,int> end_positions; // Where a given site ends; the allele (even) marker is stored.
    par_map parental_map;

private:
    std::string output_file;
    int_vec my_PRG_string;
    /**
     * Discover where site boundaries lie, and convert any odd end markers to even end markers
     */
    void map_and_normalise_ends();
};
