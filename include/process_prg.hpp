//
// Created by robyn on 09/06/17.
//

#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"

#ifndef GRAMTOOLS_PROCESS_PRG_HPP
#define GRAMTOOLS_PROCESS_PRG_HPP


using WavletTree = sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>>
using FM_Index = sdsl::csa_wt<WavletTree, 2, 16777216>;

FM_Index construct_fm_index(std::string fname,
                            std::string int_al_fname,
                            std::string memory_log_fname,
                            std::string csa_file,
                            bool fwd, bool verbose);

std::string read_prg_file(std::string prg_fpath);


#endif //GRAMTOOLS_PROCESS_PRG_HPP
