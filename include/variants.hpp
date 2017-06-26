#include "process_prg.hpp"


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


using MarkerIndex = uint64_t;
using MarkerValue = uint64_t;
using MarkerPosition = std::pair<MarkerIndex, MarkerValue>;

class MarkerPositions : public std::iterator<std::input_iterator_tag, MarkerPosition> {
public:
    MarkerPosition position;
    MarkerPosition next_position;

    MarkerPositions(unsigned long start_idx, unsigned long end_idx,
                    const FM_Index &fm, const VariantMarkers &var) : fm_index(fm), variants(var) {

        MarkerPositions::count_markers_before_start = variants.rank(start_idx);
        MarkerPositions::current_offset = 0;

        MarkerPositions::start_idx = start_idx;
        MarkerPositions::end_idx = end_idx;
    }

    MarkerPositions(const MarkerPositions &) = default;

    MarkerPositions &operator++() {
        ++current_offset;
        position = next_position;
        next_position = get_position(current_offset + 1);
        return *this;
    }

    MarkerPositions &begin() {
        current_offset = 1;
        position = get_position(current_offset);
        next_position = get_position(current_offset + 1);
        return *this;
    }

    MarkerPositions end() {
        MarkerPositions tmp(*this);
        tmp.position = std::make_pair(-1, -1);
        return tmp;
    }

    bool is_second_to_last() {
        return (next_position.first = -1) && (next_position.second = -1);
    }

    bool operator==(const MarkerPositions &rhs) const {
        return position == rhs.position;
    }

    bool operator!=(const MarkerPositions &rhs) const {
        return position != rhs.position;
    }

    MarkerPosition &operator*() {
        return position;
    }

private:
    MarkerPosition get_position(const unsigned int off) const {
        if (count_markers_before_start + off > variants.count_set_bits)
            return std::make_pair(-1, -1);

        auto marker_idx = variants.select(count_markers_before_start + off);
        if (marker_idx > end_idx)
            return std::make_pair(-1, -1);

        auto marker = fm_index.bwt[marker_idx];
        return std::make_pair(marker_idx, marker);
    }

    const FM_Index &fm_index;
    const VariantMarkers &variants;

    unsigned int count_markers_before_start;
    unsigned int current_offset;
    unsigned int start_idx, end_idx;
};


#endif //GRAMTOOLS_VARIANTS_HPP
