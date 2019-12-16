/**
 * @file
 * Code to produce the data structures supporting backward searching (FM-index using the SDSL library; bit masks) &
 * coverage recording (coverage_Graph)
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


    /**
     * Produce FM index from integer-encoded prg.
     * FM index is built using sdsl library.
     * Memory footprint of index construction is logged to disk.
     */
    FM_Index generate_fm_index(const Parameters &parameters);

    FM_Index load_fm_index(const Parameters &parameters);

    coverage_Graph generate_cov_graph(const Parameters &parameters, PRG_String const &prg_string);


    /**************
     * Bit masks **
     **************/

    /**
     * One bit vector per nucleotide in the BWT of the linearised PRG.
     * We use this to avoid rank/select queries on the BWT itself, which has an extended alphabet due to variant markers.
     */
    struct DNA_BWT_Masks {
        sdsl::bit_vector mask_a;
        sdsl::bit_vector mask_c;
        sdsl::bit_vector mask_g;
        sdsl::bit_vector mask_t;
    };

    /**
     * Generate BWT bit vector masks for each of A,C,G and T in the BWT of the prg.
     */
    DNA_BWT_Masks generate_bwt_masks(FM_Index const& fm_index,
                                     const Parameters &parameters);

    DNA_BWT_Masks load_dna_bwt_masks(const FM_Index &fm_index,
                                     const Parameters &parameters);


    /**
     * Bit vector for variant marker presence in the BWT of the prg.
     * @param fm_index which contains the bwt characters.
     */
    sdsl::bit_vector generate_bwt_markers_mask(const FM_Index &fm_index);

}

#endif //GRAMTOOLS_MK_DS_HPP
