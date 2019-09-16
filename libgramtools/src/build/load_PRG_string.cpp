#include <build/load_PRG_string.hpp>

include "load_PRG_string.hpp"

PRG_String::PRG_String(std::string const& file_in){
    fstream input(file_in, ios::in | ios::binary);
    if (!input) throw "Protobuf file not found";
    else if (!my_PRG_string.ParseFromInstream(&input)) throw "Failed to read protobuf file";
}
