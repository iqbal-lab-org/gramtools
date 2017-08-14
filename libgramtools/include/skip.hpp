#ifndef GRAMTOOLS_SKIP_HPP
#define GRAMTOOLS_SKIP_HPP


bool process_variant_edge_marker(uint64_t &left, uint64_t &right, const uint64_t maxx, const uint64_t marker_value,
                                 const FM_Index &fm_index);

bool skip(uint64_t &left, uint64_t &right, const uint64_t maxx, const uint64_t num, const FM_Index &fm_index);


#endif //GRAMTOOLS_SKIP_HPP
