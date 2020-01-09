#ifndef PRG_TYPES_HPP
#define PRG_TYPES_HPP

#include <boost/shared_ptr.hpp>

//Forward-declaration for making coverage_Graph available
class coverage_Graph;


// Forward-declarations for aliasing
class coverage_Node;
struct node_access;
struct targeted_marker;

namespace gram{
    using seqPos = int32_t;
    using covG_ptr = boost::shared_ptr<coverage_Node>;
    using marker_to_node = std::unordered_map<Marker, covG_ptr>;
    using access_vec = std::vector<node_access>;
    using target_m = std::unordered_map<Marker, std::vector<targeted_marker>>;
    using covG_ptr_map = std::map<covG_ptr, covG_ptr, std::greater<covG_ptr> >;
}

#endif //PRG_TYPES_HPP
