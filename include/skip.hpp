#include "bwt_search.h"


#ifndef GRAMTOOLS_SKIP_HPP
#define GRAMTOOLS_SKIP_HPP

bool process_variant_edge_marker(CSA &fm_index,
                                 uint64_t &left, uint64_t &right,
                                 uint64_t &left_rev, uint64_t &right_rev,
                                 const uint64_t marker_value, const uint64_t maxx);

bool skip(CSA &fm_index,
          uint64_t &left, uint64_t &right,
          uint64_t &left_rev, uint64_t &right_rev,
          uint64_t num, uint64_t maxx);

#endif //GRAMTOOLS_SKIP_HPP
