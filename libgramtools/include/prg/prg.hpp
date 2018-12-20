/**@file
 * Defines the prg-related data structure supporting quasimapping.
 * Also defines routines manipulating or populating this data structure (eg for integer encoding a prg); see gram::commands::build
 * Also defines routine for loading a prg from disk: see load_prg_info().
 */
#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>
#include <tuple>

#include "common/utils.hpp"
#include "dna_ranks.hpp"
#include "fm_index.hpp"


#ifndef GRAMTOOLS_PRG_HPP
#define GRAMTOOLS_PRG_HPP

namespace gram {

    /**
     * The key data structure holding all of the information used for vBWT backward search.
     */
    struct PRG_Info {
        FM_Index fm_index; /**< @note Accessing this data structure ([]) accesses the suffix array. */
        sdsl::int_vector<> encoded_prg;

        sdsl::int_vector<> sites_mask;
        sdsl::int_vector<> allele_mask;

        sdsl::bit_vector bwt_markers_mask; /**< Bit vector flagging variant site marker presence in bwt.*/
        sdsl::rank_support_v<1> bwt_markers_rank; /**< Used for searching for variant site markers inside an SA interval. @see `gram::left_markers_search()`.*/
        sdsl::select_support_mcl<1> bwt_markers_select; /**< Used for counting all variant site markers up to the left of an SA interval. @see `gram::left_markers_search()`.*/
        uint64_t markers_mask_count_set_bits;

        sdsl::bit_vector prg_markers_mask;
        sdsl::rank_support_v<1> prg_markers_rank;
        sdsl::select_support_mcl<1> prg_markers_select;

        DNA_BWT_Masks dna_bwt_masks; /**<Holds bit masks over the bwt for dna nucleotides. Used for rank queries to BWT during backward search. */
        // Rank support will be provided for each bit mask held in structure above.
        sdsl::rank_support_v<1> rank_bwt_a;
        sdsl::rank_support_v<1> rank_bwt_c;
        sdsl::rank_support_v<1> rank_bwt_g;
        sdsl::rank_support_v<1> rank_bwt_t;

        uint64_t max_alphabet_num;
    };

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
     * Finds largest integer in the (integer-encoded) prg.
     */
    uint64_t get_max_alphabet_num(const sdsl::int_vector<> &encoded_prg);

    /**
     * Calls prg encoding routine and stores encoded prg to file.
     * @see parse_raw_prg_file()
     */
    sdsl::int_vector<> generate_encoded_prg(const Parameters &parameters);

    sdsl::int_vector<> parse_raw_prg_file(const std::string &prg_fpath);

    /**
    * Read in file containing prg, as stream of characters.
    */
    std::string load_raw_prg(const std::string &prg_fpath);

    /**
     * Convert prg as string of characters to vector of integers.
     * Nucleotides encoded as 1-4. Variant markers can make up several characters so are treated with a buffer.
     */
    sdsl::int_vector<> encode_prg(const std::string &prg_raw);

    /**
     * Write out marker digits to the encoded prg as a single integer.
     * @see concat_marker_digits()
     */
    void flush_marker_digits(std::vector<int> &marker_digits,
                             sdsl::int_vector<> &encoded_prg,
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
