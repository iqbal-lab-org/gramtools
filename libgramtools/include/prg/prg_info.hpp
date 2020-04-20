/**@file
 * Defines a prg-related data structure holding all the structures supporting quasimapping,
 * except for the `kmer_index`.
 */
#ifndef GRAMTOOLS_PRG_INFO_HPP
#define GRAMTOOLS_PRG_INFO_HPP

#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>

#include "prg/make_data_structures.hpp"
#include "prg/linearised_prg.hpp"

namespace gram {

    /**
     * The key data structure holding all of the information used for vBWT backward search.
     */
    struct PRG_Info {
        FM_Index fm_index; /**< FM_index as a `sdsl::csa_wt` from the `sdsl` library. @note Accessing this data structure ([]) accesses the suffix array. */
        marker_vec encoded_prg;
        std::unordered_map<Marker,int> last_allele_positions;
        
        mutable coverage_Graph coverage_graph; // Can pass PRG_Info as const but still mutate this (record pb coverage)

        sdsl::bit_vector bwt_markers_mask; /**< Bit vector flagging variant site marker presence in bwt.*/
        uint64_t markers_mask_count_set_bits;

        DNA_BWT_Masks dna_bwt_masks; /**<Holds bit masks over the bwt for dna nucleotides. Used for rank queries to BWT during backward search. */
        // Rank support will be provided for each bit mask held in structure above.
        sdsl::rank_support_v<1> rank_bwt_a;
        sdsl::rank_support_v<1> rank_bwt_c;
        sdsl::rank_support_v<1> rank_bwt_g;
        sdsl::rank_support_v<1> rank_bwt_t;

        uint64_t num_variant_sites;


        // Only used for kmer indexing without `all-kmers`
        sdsl::int_vector<> sites_mask; /**< Stores the site number at each allele position. Variant markers and outside variant sites get 0.*/
        sdsl::int_vector<> allele_mask; /**< Stores the allele index at each allele position. Variant markers and outside variant sites get 0. */
        sdsl::bit_vector prg_markers_mask; /**< Bit vector flagging variant site marker presence in prg.*/
        sdsl::rank_support_v<1> prg_markers_rank;
        sdsl::select_support_mcl<1> prg_markers_select;
    };

    /**
     * Populates PRG_Info struct from disk.
     * Contains encoded prg, fm_index and masks the BWT of the prg with rank and select support.
     * Note that the fm_index contains the bwt, and that **it** has rank support.
     * @see PRG_Info()
     */
    PRG_Info load_prg_info(CommonParameters const &parameters);

}

#endif //GRAMTOOLS_PRG_INFO_HPP
