#include "variants.hpp"


VariantMarkers parse_variants(const FM_Index &fm_index){
    VariantMarkers variants;
    variants.mask = sdsl::bit_vector(fm_index.bwt.size(), 0);

    sdsl::bit_vector variant_mask(fm_index.bwt.size(), 0);
    for(unsigned int i=0; i < fm_index.bwt.size(); i++)
        variants.mask[i] = fm_index.bwt[i] > 4;

    variants.rank = sdsl::rank_support_v<1> (&variants.mask);
    variants.select = sdsl::select_support_mcl<1> (&variants.mask);
    variants.count_set_bits = variants.rank(variants.mask.size());

    return variants;
}


std::vector<std::pair<uint64_t, uint64_t>> find_variant_markers(unsigned long start_idx, unsigned long end_idx,
                                                                const FM_Index &fm_index,
                                                                const VariantMarkers &variants){

    const auto count_markers_before_start = variants.rank(start_idx);
    std::vector<std::pair<uint64_t, uint64_t>> results;

    for (int i = 1; count_markers_before_start + i <= variants.count_set_bits; ++i) {
        auto marker_idx = variants.select(count_markers_before_start + i);
        if (marker_idx > end_idx)
            break;

        auto marker = fm_index.bwt[marker_idx];
        auto result = std::make_pair(marker_idx, marker);
        results.push_back(result);
    }
    return results;
}