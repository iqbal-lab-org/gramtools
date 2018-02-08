#include <vector>
#include <string>
#include <fstream>

#include <sdsl/vectors.hpp>

#include "common/parameters.hpp"
#include "fm_index.hpp"
#include "common/utils.hpp"


#ifndef GRAMTOOLS_MASKS_H
#define GRAMTOOLS_MASKS_H

sdsl::int_vector<> load_allele_mask(const Parameters &parameters);

sdsl::int_vector<> generate_allele_mask(const sdsl::int_vector<> &encoded_prg);

sdsl::int_vector<> load_sites_mask(const Parameters &parameters);

sdsl::int_vector<> generate_sites_mask(const sdsl::int_vector<> &encoded_prg);

sdsl::bit_vector generate_prg_markers_mask(const sdsl::int_vector<> &encoded_prg);

sdsl::bit_vector generate_bwt_markers_mask(const FM_Index &fm_index);

#endif //GRAMTOOLS_MASKS_H
