#include <cstdint>
#include <string>
#include <iostream>

#include "sequence_read/seqread.hpp"
#include "common/utils.hpp"

using namespace gram;

/******************************************
 * Characters to integers and vice-versa **
 ******************************************/
EncodeResult gram::encode_char(const char &c) {
    EncodeResult encode_result = {};

    switch (c) {
        case 'A':
        case 'a':
            encode_result.is_dna = true;
            encode_result.character = 1;
            return encode_result;

        case 'C':
        case 'c':
            encode_result.is_dna = true;
            encode_result.character = 2;
            return encode_result;

        case 'G':
        case 'g':
            encode_result.is_dna = true;
            encode_result.character = 3;
            return encode_result;

        case 'T':
        case 't':
            encode_result.is_dna = true;
            encode_result.character = 4;
            return encode_result;

        default:
            /// The character is non-DNA, so must be a variant marker.
            encode_result.is_dna = false;
            encode_result.character = (uint32_t) c - '0';
            return encode_result;
    }
}

int_Base gram::encode_dna_base(const char &base_str) {
    EncodeResult result = encode_char(base_str);
    if (result.is_dna) return result.character;
    else return 0;
}

std::string  gram::decode_dna_base(const int_Base& base){
    switch(base){
        case 1: return "A";
        case 2: return "C";
        case 3: return "G";
        case 4: return "T";
        default: std::cerr << "Error: the argument " << base << " is not in [1,4]"; exit(1);
    }
}

Sequence gram::encode_dna_bases(const std::string &dna_str) {
    Sequence pattern;
    for (const auto &base_str: dna_str) {
        int_Base encoded_base = encode_dna_base(base_str);
        if (encoded_base == 0)
            return Sequence {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}


Sequence gram::encode_dna_bases(const GenomicRead &read_sequence) {
    const auto sequence_length = read_sequence.seq.size();
    Sequence pattern;
    for (uint64_t i = 0; i < sequence_length; i++) {
        int_Base encoded_base = encode_dna_base(read_sequence.seq[i]);
        if (encoded_base == 0)
            return Sequence {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}

