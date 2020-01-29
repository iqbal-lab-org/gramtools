/** @file
 * Defines a set of data structures and functions common to the backend.
 */
#ifndef GRAMTOOLS_UTILS_HPP
#define GRAMTOOLS_UTILS_HPP

#include <vector>
#include <list>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>
#include "sequence_read/seqread.hpp"


namespace gram {

    std::string full_path(const std::string &gram_dirpath,
                          const std::string &file_name);

    /******************
     * Data typedefs **
     ******************/
    using int_Base = uint8_t; /**< nucleotide represented as byte-sized integer */
    using Sequence = std::vector<int_Base>; /** A string of nucleotides is represented as a vector of `Base`s. */
    using Sequences = std::vector<Sequence>;

    using Marker = uint32_t ; /**< An integer >=5 describing a site or allele marker in the prg. */
    using marker_vec = std::vector<Marker>;
    using AlleleId = uint32_t; /**< An integer describing which allele is referred to within a given variant site. */
    using AlleleIds = std::vector<AlleleId>;
    using AlleleIdSet = std::set<AlleleId>;
    using VariantLocus = std::pair<Marker, AlleleId>; /**< A Variant site/`AlleleId` combination.*/

    using parental_map = std::unordered_map<Marker, VariantLocus>; /** Map of a site to its parental Locus */
    using haplo_map = std::unordered_map<AlleleId, marker_vec>;
    /**
     * Opposite of parental_map. Used for associating a site's haplogroup to IDs of sites sitting inside it.
     * The top-level key is a site ID, and you then access each haplogroup in the inner map.
     */
    using child_map = std::unordered_map<Marker, haplo_map >;

    // coverage-related
    using CovCount = uint16_t;
    using PerBaseCoverage = std::vector<CovCount>; /**< Number of reads mapped to each base of an allele */
    using PerAlleleCoverage = std::vector<CovCount>; /**< Number of reads mapped to each of several alleles */

    /**
     * Build child_map from parental_map
     */
     child_map build_child_map(parental_map const& par_map);

    /**
     * Unified definition of a site marker
     */
    bool is_site_marker(Marker const& variant_marker);

    /**
     * Unified definition of an allele marker
     */
    bool is_allele_marker(Marker const& variant_marker);

    void ensure_is_site_marker(Marker const& site_ID);

    /**
     * Conversion of a site ID to a 0-based index, suitable for site array access.
     * In other words maps 5 to 0, 7 to 1, etc.
     */
     std::size_t siteID_to_index(Marker const& site_ID);

    /**
     * Opposite conversion: 0-based array access index to site ID
     */
    Marker index_to_siteID(std::size_t const& idx);

    /******************************************
     * Characters to integers and vice-versa **
     ******************************************/
    struct EncodeResult {
        bool is_dna;
        uint32_t character;
    };

    /**
     * Encode a single character read from the linearised prg as an integer.
     * Use `EncodeResult` object to additionally store if the encoded character is DNA or not.
     * @see EncodeResult()
     */
    EncodeResult encode_char(const char &c);

    /**
     * Encode dna character as integer (range: 1-4)
     */
    int_Base encode_dna_base(const char &base_str);

    /**
     * Decode integer into dna base, as std::string.
     */
    std::string decode_dna_base(const int_Base& base);

    Sequence encode_dna_bases(const std::string &dna_str);

    Sequence encode_dna_bases(const GenomicRead &read_sequence);


    /************************************************
     * Linearised PRGs to integer vectors, and v-v.**
     ************************************************/

    /** Converts linearised PRG as int vector to a more readable string.
     *  We use the following notation: '[' opens a site, ',' delimits alleles in a site, ']' closes a site.
     * */
    std::string ints_to_prg_string(std::vector<Marker> const& int_vec);


    /**
     * Convert a nested PRG string to int representation, with linear site numbering.
     * The site numbering is based on the fixed order in which '[' chars are encountered;
     * thus int -> prg_string -> int can lose original site numbering.
     */
    std::vector<Marker> prg_string_to_ints(std::string const& string_prg);


    /************
     * Hashing **
     ************/
    template<typename SEQUENCE>
    struct sequence_hash {
        std::size_t operator()(const SEQUENCE &seq) const {
            std::size_t hash = 0; // Used as an initial seed
            boost::hash_range(hash, seq.begin(), seq.end());
            return hash;
        }
    };

    template<typename SEQUENCE, typename T>
    using SequenceHashMap = std::unordered_map<SEQUENCE, T, sequence_hash<SEQUENCE>>;

    template<typename PAIR, typename T>
    using PairHashMap = std::unordered_map<PAIR, T, boost::hash<PAIR>>;

    template<typename T>
    using HashSet = std::unordered_set<T, boost::hash<T>>;

}

#endif //GRAMTOOLS_UTILS_HPP
