#include <cstdint>
#include <vector>
#include <string>
#include <iostream>

#include "utils.hpp"


uint8_t encode_dna_base(const char &base_str) {
    switch (base_str) {
        case 'A':
        case 'a':
            return 1;

        case 'C':
        case 'c':
            return 2;

        case 'G':
        case 'g':
            return 3;

        case 'T':
        case 't':
            return 4;

        default:
            std::cout << "Error encoding base" << std::endl;
            break;
    }
}


std::vector<uint8_t> encode_dna_bases(const std::string &dna_str) {
    std::vector<uint8_t> dna;
    for (const auto &base_str: dna_str) {
        int encoded_base = encode_dna_base(base_str);
        dna.emplace_back(encoded_base);
    }
    return dna;
}