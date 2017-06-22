#include <sdsl/suffix_arrays.hpp>
#include "bwt_search.h"

using namespace sdsl;

// TODO: rename to get_variant_site_edge?
std::pair<uint32_t, std::vector<int>> get_location(const FM_Index &fm_index,
                                                   const uint64_t marker_idx, const uint64_t marker,
                                                   const bool last, std::vector<int> &allele,
                                                   const std::vector<int> &mask_a) {
    uint64_t site_edge_marker;

    bool marker_is_site_edge = marker % 2 == 1;
    if (marker_is_site_edge) {
        site_edge_marker = marker;
        if (!last)
            allele.push_back(1);
    } else {
        site_edge_marker = marker - 1;
        allele.push_back(mask_a[fm_index[marker_idx]]);
    }
    return std::make_pair(site_edge_marker, allele);
}

//removed left and right args because they only form a proper interval is when last is true; but then they can't help with location because allele is unknown and we're just addin site

//problems:
//avoid mask_s
//when have matches to the right of both odd numbers
