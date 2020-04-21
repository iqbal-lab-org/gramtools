#ifndef GRAMTOOLS_DATA_TYPES_HPP
#define GRAMTOOLS_DATA_TYPES_HPP

#include <cstdint>
#include <vector>
#include <set>

#define FIRST_ALLELE 0
#define ALLELE_UNKNOWN -1 // This signifier must NEVER be a possible allele ID

namespace gram {
    using int_Base = uint8_t; /**< nucleotide represented as byte-sized integer */
    using Sequence = std::vector<int_Base>; /** A string of nucleotides is represented as a vector of `Base`s. */
    using Sequences = std::vector<Sequence>;

    using Marker = uint32_t; /**< An integer >=5 describing a site or allele marker in the prg. */
    using marker_vec = std::vector<Marker>;
    using AlleleId = int32_t; /**< An integer describing which allele is referred to within a given variant site. */
    using AlleleIds = std::vector<AlleleId>;
    using AlleleIdSet = std::set<AlleleId>;
    using VariantLocus = std::pair<Marker, AlleleId>; /**< A Variant site/`AlleleId` combination.*/

    static bool is_site_marker(Marker const& variant_marker){
        if (!(variant_marker > 4)) throw std::invalid_argument("The given marker is not a variant marker (>4)");
        return variant_marker % 2 == 1;
    }

    static bool is_allele_marker(Marker const& variant_marker){
        return !is_site_marker(variant_marker);
    }

    static void ensure_is_site_marker(Marker const& site_ID){
        if (!is_site_marker(site_ID)) throw std::invalid_argument("The given marker is not a site ID");
    }

    /**
     * Conversion of a site ID to a 0-based index, suitable for site array access.
     * In other words maps 5 to 0, 7 to 1, etc.
     */
    static std::size_t siteID_to_index(Marker const& site_ID){
        ensure_is_site_marker(site_ID);
        return (site_ID - 5) / 2;
    }

    /**
     * Opposite conversion: 0-based array access index to site ID
     */
    static Marker index_to_siteID(std::size_t const& idx){
        return idx * 2 + 5;
    }
}

#endif //GRAMTOOLS_DATA_TYPES_HPP
