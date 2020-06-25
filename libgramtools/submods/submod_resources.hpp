#ifndef COMMON_GEN_PRG
#define COMMON_GEN_PRG

#include "prg/prg_info.hpp"

namespace gram::submods {
gram::PRG_Info generate_prg_info(const marker_vec &prg_raw);

std::string decode(uint64_t base);

using covG_ptrPair = std::pair<covG_ptr, covG_ptr>;

/**
 * Given a map of all bubbles and a `siteID` of interest, returns the pair of
 * `covG_ptr` corresponding to the start and end nodes of the site.
 */
covG_ptrPair get_bubble_nodes(covG_ptr_map bubble_map, Marker site_ID);
}  // namespace gram::submods

#endif  // COMMON_GEN_PRG
