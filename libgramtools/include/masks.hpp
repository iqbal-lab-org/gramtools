#include <vector>
#include <string>
#include <fstream>

#include <sdsl/vectors.hpp>

#include "parameters.hpp"
#include "fm_index.hpp"
#include "utils.hpp"


#ifndef GRAMTOOLS_MASKS_H
#define GRAMTOOLS_MASKS_H

sdsl::int_vector<> load_allele_mask(const Parameters &parameters);

sdsl::int_vector<> generate_allele_mask(const sdsl::int_vector<> &encoded_prg);

sdsl::bit_vector generate_prg_markers_mask(const sdsl::int_vector<> &encoded_prg);

sdsl::bit_vector generate_bwt_markers_mask(const FM_Index &fm_index);

class MasksParser {
public:
    uint64_t max_alphabet_num;
    std::vector<uint64_t> sites;
    MasksParser() {};
    MasksParser(const std::string &sites_fname);
private:
    void parse_sites(std::istream &stream);
};

#endif //GRAMTOOLS_MASKS_H
