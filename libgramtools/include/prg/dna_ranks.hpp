/**
 * @file
 * Functions for producing, storing and loading bit masks over the BWT of the prg.
 * A bit mask is a bit vector with a bit set if a given nucleotide base is present in the BWT, and unset otherwise.
 */
#include "fm_index.hpp"


#ifndef GRAMTOOLS_DNA_RANKS_HPP
#define GRAMTOOLS_DNA_RANKS_HPP

namespace gram {

    struct DNA_BWT_Masks {
        sdsl::bit_vector mask_a;
        sdsl::bit_vector mask_c;
        sdsl::bit_vector mask_g;
        sdsl::bit_vector mask_t;
    };

    /**
     * Generate BWT bit vector masks for each of A,C,G and T in the BWT of the prg.
     */
    void generate_dna_bwt_masks(const FM_Index &fm_index,
                                const Parameters &parameters);

    DNA_BWT_Masks load_dna_bwt_masks(const FM_Index &fm_index,
                                     const Parameters &parameters);

}

#endif //GRAMTOOLS_DNA_RANKS_HPP
