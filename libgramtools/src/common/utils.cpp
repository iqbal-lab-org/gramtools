#include <cstdint>
#include <vector>
#include <string>
#include <iostream>

#include <boost/filesystem.hpp>

#include "sequence_read/seqread.hpp"
#include "common/utils.hpp"


namespace fs = boost::filesystem;


std::string full_path(const std::string &gram_dirpath,
                      const std::string &file_name) {
    fs::path dir(gram_dirpath);
    fs::path file(file_name);
    fs::path full_path = dir / file;
    return full_path.string();
}


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


Base compliment_encoded_base(const Base &encoded_base) {
    switch (encoded_base) {
        case 1:
            return 4;
        case 2:
            return 3;
        case 3:
            return 2;
        case 4:
            return 1;
        default:
            return 0;
    }
}


Pattern reverse_compliment_read(const Pattern &read) {
    Pattern reverse_read;
    reverse_read.reserve(read.size());

    for (auto it = read.rbegin(); it != read.rend(); ++it) {
        const auto &base = *it;
        auto compliment_base = compliment_encoded_base(base);
        reverse_read.push_back(compliment_base);
    }
    return reverse_read;
}


Pattern encode_dna_bases(const GenomicRead &read_sequence) {
    const auto sequence_length = strlen(read_sequence.seq);
    Pattern pattern;
    for (uint64_t i = 0; i < sequence_length; i++) {
        Base encoded_base = encode_dna_base(read_sequence.seq[i]);
        if (encoded_base == 0)
            return Pattern {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}