/**
 * @file
 * Generates an FM-index of an encoded prg using the SDSL library.
 */
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "common/parameters.hpp"


#ifndef GRAMTOOLS_PROCESS_PRG_HPP
#define GRAMTOOLS_PROCESS_PRG_HPP


namespace gram {

    using WaveletTree = sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>>;
    using FM_Index = sdsl::csa_wt<WaveletTree, 1, 16777216>; /**< The two numbers are the sampling densities for SA and ISA. 1 means all SA entries are stored.*/


    FM_Index load_fm_index(const Parameters &parameters);

    /**
     * Produce FM index from integer-encoded prg.
     * FM index is built using sdsl library.
     * Memory footprint of index construction is logged to disk.
     */
    FM_Index generate_fm_index(const Parameters &parameters);

}

#endif //GRAMTOOLS_PROCESS_PRG_HPP
