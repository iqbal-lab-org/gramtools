/**
 * @file
 * Code to produce the data structures supporting backward searching (FM-index
 * using the SDSL library; bit masks) & coverage recording (coverage_Graph)
 */

#ifndef GRAMTOOLS_MK_DS_HPP
#define GRAMTOOLS_MK_DS_HPP

#include "build/parameters.hpp"
#include "prg/linearised_prg.hpp"
#include "prg/types.hpp"

namespace gram {
/**
 * Produce FM index from integer-encoded prg.
 * FM index is built using sdsl library.
 * Memory footprint of index construction is logged to disk.
 */
FM_Index generate_fm_index(BuildParams const &parameters);

FM_Index load_fm_index(CommonParameters const &parameters);

/**************
 * Cov graph***
 **************/

coverage_Graph generate_cov_graph(CommonParameters const &parameters,
                                  PRG_String const &prg_string);

/**
 * Build child_map from parental_map
 */
child_map build_child_map(parental_map const &par_map);

/**************
 * Bit masks **
 **************/

/**
 * Generate BWT bit vector masks for each of A,C,G and T in the BWT of the prg.
 */
DNA_BWT_Masks generate_bwt_masks(FM_Index const &fm_index,
                                 CommonParameters const &parameters);

DNA_BWT_Masks load_dna_bwt_masks(const FM_Index &fm_index,
                                 CommonParameters const &parameters);

/**
 * Bit vector for variant marker presence in the BWT of the prg.
 * @param fm_index which contains the bwt characters.
 */
sdsl::bit_vector generate_bwt_markers_mask(const FM_Index &fm_index);

}  // namespace gram

#endif  // GRAMTOOLS_MK_DS_HPP
