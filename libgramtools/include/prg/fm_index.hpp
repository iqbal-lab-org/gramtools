#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "common/parameters.hpp"


#ifndef GRAMTOOLS_PROCESS_PRG_HPP
#define GRAMTOOLS_PROCESS_PRG_HPP

using WavletTree = sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>>;
using FM_Index = sdsl::csa_wt<WavletTree, 2, 16777216>;

// skip optimisation (need to rebuild fm-index):
// using FM_Index = sdsl::csa_wt<WavletTree, 2, 2>;

FM_Index load_fm_index(const Parameters &parameters);

FM_Index generate_fm_index(const Parameters &parameters);

#endif //GRAMTOOLS_PROCESS_PRG_HPP
