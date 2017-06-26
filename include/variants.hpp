#include "fm_index.hpp"


#ifndef GRAMTOOLS_VARIANTS_HPP
#define GRAMTOOLS_VARIANTS_HPP


struct VariantMarkers {
    sdsl::bit_vector mask;
    sdsl::rank_support_v<1> rank;
    sdsl::select_support_mcl<1> select;
    unsigned long count_set_bits;
};

VariantMarkers parse_variants(const FM_Index &fm_index);

std::vector<std::pair<uint64_t, uint64_t>> find_variant_markers(unsigned long start_idx, unsigned long end_idx,
                                                                const FM_Index &fm_index,
                                                                const VariantMarkers &variant_mask);

#endif //GRAMTOOLS_VARIANTS_HPP
