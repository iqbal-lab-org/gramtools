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

using MarkerIndex = uint64_t;
using MarkerValue = uint64_t;
using MarkerPosition = std::pair<MarkerIndex, MarkerValue>;

class MarkerPositions : public std::iterator<std::input_iterator_tag, MarkerPosition> {
public:
    MarkerPosition position;
    MarkerPosition next_position;

    MarkerPositions(unsigned long start_idx, unsigned long end_idx,
                    const FM_Index &fm, const VariantMarkers &var);
    MarkerPositions(const MarkerPositions &) = default;
    MarkerPositions &operator++();
    MarkerPositions &begin();
    MarkerPositions end() const;
    bool is_second_to_last() const;

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
    MarkerPosition get_position(const unsigned int off) const;
    const FM_Index &fm_index;
    const VariantMarkers &variants;

    unsigned int count_markers_before_start;
    unsigned int current_offset;
    unsigned int start_idx, end_idx;
};


#endif //GRAMTOOLS_VARIANTS_HPP
