#ifndef PRG_TYPES_HPP
#define PRG_TYPES_HPP

#include <map>

#include <boost/shared_ptr.hpp>

#include "common/data_types.hpp"

// Forward-declaration for making coverage_Graph available
class coverage_Graph;

// Forward-declarations for aliasing
class coverage_Node;
struct node_access;
struct targeted_marker;

namespace gram {
using covG_ptr = boost::shared_ptr<coverage_Node>;
using marker_to_node = std::unordered_map<Marker, covG_ptr>;
using access_vec = std::vector<node_access>;
using target_m = std::unordered_map<Marker, std::vector<targeted_marker>>;
using covG_ptr_map = std::map<covG_ptr, covG_ptr, std::greater<covG_ptr>>;

using parental_map =
    std::unordered_map<Marker,
                       VariantLocus>; /** Map of a site to its parental Locus */
using haplo_map = std::unordered_map<AlleleId, marker_vec>;
/**
 * Opposite of parental_map. Used for associating a site's haplogroup to IDs of
 * sites sitting inside it. The top-level key is a site ID, and you then access
 * each haplogroup in the inner map.
 */
using child_map = std::unordered_map<Marker, haplo_map>;
}  // namespace gram

#endif  // PRG_TYPES_HPP
