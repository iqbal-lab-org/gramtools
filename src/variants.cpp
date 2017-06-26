#include "variants.hpp"


VariantMarkers parse_variants(const FM_Index &fm_index) {
    VariantMarkers variants;
    variants.mask = sdsl::bit_vector(fm_index.bwt.size(), 0);

    sdsl::bit_vector variant_mask(fm_index.bwt.size(), 0);
    for (unsigned int i = 0; i < fm_index.bwt.size(); i++)
        variants.mask[i] = fm_index.bwt[i] > 4;

    variants.rank = sdsl::rank_support_v<1>(&variants.mask);
    variants.select = sdsl::select_support_mcl<1>(&variants.mask);
    variants.count_set_bits = variants.rank(variants.mask.size());

    return variants;
}


std::vector<std::pair<uint64_t, uint64_t>> find_variant_markers(unsigned long start_idx, unsigned long end_idx,
                                                                const FM_Index &fm_index,
                                                                const VariantMarkers &variants) {

    MarkerPositions marker_positions(start_idx, end_idx, fm_index, variants);

    std::vector<std::pair<uint64_t, uint64_t>> results;
    for (auto marker_position = marker_positions.begin();
         marker_position != marker_positions.end();
         ++marker_position) {

        results.push_back(marker_position.position);
    }
    return results;
}