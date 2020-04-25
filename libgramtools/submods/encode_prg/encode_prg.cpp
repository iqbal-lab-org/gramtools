/**
 * @file Takes a linearised prg encoded as follows:
 *  - '[' and ']' mark site start and site end points
 *  - ',' separates alleles in a site
 *  - {A,C,G,T} for normal sequence
 *
 *  And turns it into a vector of integers representing the prg,
 *  suitable as --prg parameter to `gramtools build` command.
 */

#include <iostream>

#include "prg/linearised_prg.hpp"

void usage(const char* argv[]){
    std::cout << "Usage: " << argv[0] << " -o PATH" << std::endl;
    std::cout << "Provide input sequence on stdin, can be fasta" << std::endl;
    std::cout << "The prg string variant markers should be encoded using '[', ']' and ','" << std::endl;
    std::cout << "View output using:  hexdump -v -e '1/4 \"%d \"' " << std::endl;
    exit(1);
}

int main(int argc, const char* argv[]) {
    if (argc != 3) usage(argv);

    std::ofstream fout;
    std::string flag{argv[1]};
    if (flag != "-o") usage(argv);
    std::string fout_path = argv[2];

    std::string line, prg_string;
    while (std::cin){
        std::getline(std::cin, line);
        if (line[0] == '>') continue;
        std::transform(line.begin(), line.end(), line.begin(), ::toupper);
        prg_string += line;
    }

    if (prg_string.empty()) usage(argv);

    auto as_marker_vec = prg_string_to_ints(prg_string);

    PRG_String handler(as_marker_vec);
    handler.write(fout_path);

    std::cout << "Made integer-encoded linearised PRG." << std::endl;
    std::cout << "Use:  hexdump -v -e '1/4 \"%d \"' " << fout_path <<
              " to get textual representation." << std::endl;
}
