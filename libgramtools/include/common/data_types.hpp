#ifndef GRAMTOOLS_DATA_TYPES_HPP
#define GRAMTOOLS_DATA_TYPES_HPP

#include <cstdint>
#include <set>
#include <vector>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#define FIRST_ALLELE 0
#define ALLELE_UNKNOWN -1  // This signifier must NEVER be a possible allele ID

namespace gram {

using int_Base = uint8_t; /**< nucleotide represented as byte-sized integer */
using Sequence =
    std::vector<int_Base>; /** A string of nucleotides is represented as a
                              vector of `Base`s. */
using Sequences = std::vector<Sequence>;

using Marker = uint32_t; /**< An integer >=5 describing a site or allele marker
                            in the prg. */
using marker_vec = std::vector<Marker>;
using AlleleId = int32_t; /**< An integer describing which allele is referred to
                             within a given variant site. */
using AlleleIds = std::vector<AlleleId>;
using AlleleIdSet = std::set<AlleleId>;
using VariantLocus =
    std::pair<Marker, AlleleId>; /**< A Variant site/`AlleleId` combination.*/

// BWT-related
using WaveletTree = sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>>;
using FM_Index =
    sdsl::csa_wt<WaveletTree, 1,
                 16777216>; /**< The two numbers are the sampling densities for
                               SA and ISA. 1 means all SA entries are stored.*/

/**
 * One bit vector per nucleotide in the BWT of the linearised PRG.
 * We use this to avoid rank/select queries on the BWT itself, which has an
 * extended alphabet due to variant markers.
 */
struct DNA_BWT_Masks {
  sdsl::bit_vector mask_a;
  sdsl::bit_vector mask_c;
  sdsl::bit_vector mask_g;
  sdsl::bit_vector mask_t;
};

// coverage-related
using CovCount = uint16_t;
using PerBaseCoverage = std::vector<CovCount>;   /**< Number of reads mapped to
                                                    each base of an allele */
using PerAlleleCoverage = std::vector<CovCount>; /**< Number of reads mapped to
                                                    each of several alleles */

static bool is_site_marker(Marker const& variant_marker) {
  if (!(variant_marker > 4))
    throw std::invalid_argument(
        "The given marker is not a variant marker (>4)");
  return variant_marker % 2 == 1;
}

static bool is_allele_marker(Marker const& variant_marker) {
  return !is_site_marker(variant_marker);
}

static void ensure_is_site_marker(Marker const& site_ID) {
  if (!is_site_marker(site_ID))
    throw std::invalid_argument("The given marker is not a site ID");
}

/**
 * Conversion of a site ID to a 0-based index, suitable for site array access.
 * In other words maps 5 to 0, 7 to 1, etc.
 */
static std::size_t siteID_to_index(Marker const& site_ID) {
  ensure_is_site_marker(site_ID);
  return (site_ID - 5) / 2;
}

/**
 * Opposite conversion: 0-based array access index to site ID
 */
static Marker index_to_siteID(std::size_t const& idx) { return idx * 2 + 5; }
}  // namespace gram

#endif  // GRAMTOOLS_DATA_TYPES_HPP
