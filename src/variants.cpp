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

    std::vector<std::pair<uint64_t, uint64_t>> results;
    unsigned long count_pre_range = variants.rank(start_idx);
    unsigned long i = 1;

    while (true){
        if (count_pre_range + i > variants.count_set_bits)
            break;

        unsigned long idx = variants.select(count_pre_range + i);
        if (idx > end_idx)
            break;

        const auto result = std::make_pair(idx, fm_index.bwt[idx]);
        results.push_back(result);
        i++;
    }

    return results;
}