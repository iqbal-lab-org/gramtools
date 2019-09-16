#include "int_vector.pb.h"

using std::unordered_map<int, std::pair<int,int>> = par_map;

class PRG_String{
public:
    PRG_String(){};
    PRG_String(std::string const& file_in);
    PRG_String(Proto::PRG_String const& proto_in);

    void parse(); // Traverse the PRG string, figuring out the write config & mapping site parental relationships
    void write();

    /**
     * Variables
     */
     bool odd_site_end_found;
     std::set<int> start_positions;
    std::set<int> end_positions;
    par_map parental_map;
private:
    Proto::PRG_String my_PRG_string;
};
