#include <vector>
#include <string>
#include <fstream>

#include <sdsl/vectors.hpp>

#include "utils.hpp"


#ifndef GRAMTOOLS_MASKS_H
#define GRAMTOOLS_MASKS_H

sdsl::bit_vector generate_markers_mask(const sdsl::int_vector<> &encoded_prg);

class MasksParser {
public:
    uint64_t max_alphabet_num;
    std::vector<uint64_t> sites;

    std::vector<AlleleId> allele;

    MasksParser(){};
    MasksParser(const std::string &sites_fname, const std::string &alleles_fname);

private:
    void parse_sites(std::istream &stream);
    void parse_allele(std::istream &stream);
};

#endif //GRAMTOOLS_MASKS_H
