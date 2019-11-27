/**@file
 * Defines the prg-related data structure supporting quasimapping.
 * Also defines routines manipulating or populating this data structure (eg for integer encoding a prg); see gram::commands::build
 * Also defines routine for loading a prg from disk: see load_prg_info().
 */
#ifndef GRAMTOOLS_PRG_HPP
#define GRAMTOOLS_PRG_HPP

#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>
#include <tuple>

#include "common/utils.hpp"
#include "dna_ranks.hpp"
#include "prg/make_data_structures.hpp"
#include "prg/load_PRG_string.hpp"
#include "prg/coverage_graph.hpp"

namespace gram {

    using marker_map = std::unordered_map<Marker,int>;
    /**
     * The key data structure holding all of the information used for vBWT backward search.
     */
    struct PRG_Info {
        FM_Index fm_index; /**< FM_index as a `sdsl::csa_wt` from the `sdsl` library. @note Accessing this data structure ([]) accesses the suffix array. */
        marker_vec encoded_prg;
        marker_map last_allele_positions;
        
        mutable coverage_Graph coverage_graph;

        sdsl::int_vector<> sites_mask; /**< Stores the site number at each allele position. Variant markers and outside variant sites get 0.*/
        sdsl::int_vector<> allele_mask; /**< Stores the allele index at each allele position. Variant markers and outside variant sites get 0. */

        sdsl::bit_vector bwt_markers_mask; /**< Bit vector flagging variant site marker presence in bwt.*/
        uint64_t markers_mask_count_set_bits;

        sdsl::bit_vector prg_markers_mask; /**< Bit vector flagging variant site marker presence in prg.*/
        sdsl::rank_support_v<1> prg_markers_rank;
        sdsl::select_support_mcl<1> prg_markers_select;

        DNA_BWT_Masks dna_bwt_masks; /**<Holds bit masks over the bwt for dna nucleotides. Used for rank queries to BWT during backward search. */
        // Rank support will be provided for each bit mask held in structure above.
        sdsl::rank_support_v<1> rank_bwt_a;
        sdsl::rank_support_v<1> rank_bwt_c;
        sdsl::rank_support_v<1> rank_bwt_g;
        sdsl::rank_support_v<1> rank_bwt_t;

        uint64_t num_variant_sites;
    };

    /**
     * Maps the end positions of sites, producing a map with keys being the even `Marker` for that site.
     */
    marker_map map_site_ends(sdsl::int_vector<> const& encoded_prg);

    /**
     * Performs a rank query on the BWT. given a `dna_base` and a `upper_index`.
     * @param upper_index the index into the suffix array/BWT.
     * @param dna_base the base to count in the BWT.
     * @return the number of occurrences of `dna_base` up to (and excluding) `upper_index` in the BWT of the prg.
     */
    uint64_t dna_bwt_rank(const uint64_t &upper_index,
                          const Marker &dna_base,
                          const PRG_Info &prg_info);

    /**
     * Calls prg encoding routine and stores encoded prg to file.
     * @see parse_raw_prg_file()
     */
    marker_vec generate_encoded_prg(const Parameters &parameters);

    marker_vec parse_raw_prg_file(const std::string &prg_fpath);

    /**
    * Read in file containing prg, as stream of characters.
    */
    std::string load_raw_prg(const std::string &prg_fpath);

    /**
     * Convert prg as string of characters to vector of integers.
     * Nucleotides encoded as 1-4. Variant markers can make up several characters so are treated with a buffer.
     */
    marker_vec encode_prg(const std::string &prg_raw);

    /**
     * Write out marker digits to the encoded prg as a single integer.
     * @see concat_marker_digits()
     */
    void flush_marker_digits(std::vector<int> &marker_digits,
                             marker_vec &encoded_prg,
                             uint64_t &count_chars);

    /**
     * Converts a sequence of digits (0-9) into a single integer.
     */
    uint64_t concat_marker_digits(const std::vector<int> &marker_digits);

    struct EncodeResult {
        bool is_dna;
        uint32_t character;
    };

    /**
     * Encode a character read from the prg as an integer.
     * Use `EncodeResult` object to additionally store if the encoded character is DNA or not.
     * @see EncodeResult()
     */
    EncodeResult encode_char(const char &c);

    /**
     * Populates PRG_Info struct from disk.
     * Contains encoded prg, fm_index and masks over the prg and the BWT of the prg with rank and select support.
     * Note that the fm_index contains the bwt, and that **it** has rank support.
     * @see PRG_Info()
     */
    PRG_Info load_prg_info(const Parameters &parameters);

}

#endif //GRAMTOOLS_PRG_HPP
