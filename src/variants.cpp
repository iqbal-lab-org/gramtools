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


MarkerPositions::MarkerPositions(unsigned long start_idx, unsigned long end_idx,
                                 const FM_Index &fm, const VariantMarkers &var) : fm_index(fm), variants(var) {

    MarkerPositions::count_markers_before_start = variants.rank(start_idx);
    MarkerPositions::current_offset = 0;

    MarkerPositions::start_idx = start_idx;
    MarkerPositions::end_idx = end_idx;
}

MarkerPositions &MarkerPositions::operator++() {
    ++current_offset;
    position = next_position;
    next_position = get_position(current_offset + 1);
    return *this;
}

MarkerPositions &MarkerPositions::begin() {
    current_offset = 1;
    position = get_position(current_offset);
    next_position = get_position(current_offset + 1);
    return *this;
}

MarkerPositions MarkerPositions::end() const {
    MarkerPositions tmp(*this);
    tmp.position = std::make_pair(-1, -1);
    return tmp;
}

bool MarkerPositions::is_second_to_last() const {
    return (next_position.first == -1) && (next_position.second == -1);
}

MarkerPosition MarkerPositions::get_position(const unsigned int off) const {
    if (count_markers_before_start + off > variants.count_set_bits)
        return std::make_pair(-1, -1);

    auto marker_idx = variants.select(count_markers_before_start + off);
    if (marker_idx > end_idx)
        return std::make_pair(-1, -1);

    auto marker = fm_index.bwt[marker_idx];
    return std::make_pair(marker_idx, marker);
}
