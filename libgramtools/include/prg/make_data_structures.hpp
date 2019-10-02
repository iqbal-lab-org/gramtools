/**
 * @file
 * Generates an FM-index of an encoded prg using the SDSL library.
 */

#ifndef GRAMTOOLS_MK_DS_HPP
#define GRAMTOOLS_MK_DS_HPP

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <common/utils.hpp>
#include "common/parameters.hpp"
#include "coverage_graph.hpp"


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

    coverage_Graph generate_cov_graph(const Parameters &parameters, PRG_String const &prg_string);
}

#endif //GRAMTOOLS_MK_DS_HPP
