/**
 * @file
 * Routines for generating masks into the prg and the BWT of the prg.
 * The masks are vectors marking presence/absence of variant markers and site and allele numbers for nucleotides inside variant sites.
 */
#include <vector>
#include <string>
#include <fstream>

#include <sdsl/vectors.hpp>

#include "common/parameters.hpp"
#include "fm_index.hpp"
#include "common/utils.hpp"


#ifndef GRAMTOOLS_MASKS_H
#define GRAMTOOLS_MASKS_H

namespace gram {

    sdsl::int_vector<> load_allele_mask(const Parameters &parameters);

    /**
     * Generates an integer vector into the prg labelling indices within variant sites with the allele id.
     * Indices outside variant sites, and variant markers (alleles and sites) are marked as 0.
     */
    sdsl::int_vector<> generate_allele_mask(const marker_vec &encoded_prg);

    sdsl::int_vector<> load_sites_mask(const Parameters &parameters);

    /**
     * Generates an integer vector into the prg of the site number at all indices within a site.
     * At variant markers (**both** allele and site) and outside variant sites, stores a 0.
     */
    sdsl::int_vector<> generate_sites_mask(const marker_vec &encoded_prg);

    /**
     * Bit vector for variant marker presence in the prg.
     * Variant marker is a site marker (odd integer) or allele marker (even integer).
     */
    sdsl::bit_vector generate_prg_markers_mask(const marker_vec &encoded_prg);

    /**
     * Bit vector for variant marker presence in the BWT of the prg.
     * @param fm_index which contains the bwt characters.
     */
    sdsl::bit_vector generate_bwt_markers_mask(const FM_Index &fm_index);

}

#endif //GRAMTOOLS_MASKS_H
