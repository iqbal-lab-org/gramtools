#include <cstdint>
#include <vector>
#include <string>
#include <iostream>

#include "sequence_read/seqread.hpp"
#include "utils.hpp"


Base encode_dna_base(const char &base_str) {
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
            return 0;
    }
}


Pattern encode_dna_bases(const std::string &dna_str) {
    Pattern pattern;
    for (const auto &base_str: dna_str) {
        Base encoded_base = encode_dna_base(base_str);
        if (encoded_base == 0)
            return Pattern {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}


Pattern encode_dna_bases(const GenomicRead &read_sequence) {
    const auto sequence_length = strlen(read_sequence.seq);
    Pattern pattern;
    for (uint32_t i = 0; i < sequence_length; i++) {
        Base encoded_base = encode_dna_base(read_sequence.seq[i]);
        if (encoded_base == 0)
            return Pattern {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}