#include <sdsl/suffix_arrays.hpp>

#include "process_prg.hpp"
#include "skip.hpp"


bool skip(const FM_Index &fm_index,
          uint64_t &left, uint64_t &right,
          uint64_t &left_rev, uint64_t &right_rev,
          const uint64_t marker_value, const uint64_t maxx) {

    const auto is_edge_marker = marker_value % 2 == 1;
    if (is_edge_marker) {
        bool last = process_variant_edge_marker(fm_index,
                                                left, right,
                                                left_rev, right_rev,
                                                marker_value, maxx);
        return last;
    }

    const uint64_t num_begin = fm_index.C[fm_index.char2comp[marker_value - 1]];
    const auto site_boundary_lower = fm_index[num_begin];
    const auto site_boundary_upper = fm_index[num_begin + 1];

    if (site_boundary_lower < site_boundary_upper) {
        left = num_begin;
        right = num_begin + 1;

        left_rev = num_begin + 1;
        right_rev = num_begin + 2;
    } else {
        left = num_begin + 1;
        right = num_begin + 2;

        left_rev = num_begin;
        right_rev = num_begin + 1;
    }

    return false;
}

bool process_variant_edge_marker(const FM_Index &fm_index,
                                 uint64_t &left, uint64_t &right,
                                 uint64_t &left_rev, uint64_t &right_rev,
                                 const uint64_t marker_value, const uint64_t maxx) {

    const uint64_t num_begin = fm_index.C[fm_index.char2comp[marker_value]];
    const auto site_boundary_lower = fm_index[num_begin];
    const auto site_boundary_upper = fm_index[num_begin + 1];

    uint64_t ind_start;
    if (site_boundary_lower > site_boundary_upper)
        ind_start = num_begin + 1;
    else
        ind_start = num_begin;

    bool last = false;
    if (right - left == 1) {
        const uint64_t site_end = std::max(site_boundary_lower, site_boundary_upper);
        if (fm_index[left] == site_end + 1) {
            last = true;

            left = num_begin;
            if (marker_value + 2 <= maxx)
                right = fm_index.C[fm_index.char2comp[marker_value + 2]];
            else
                right = fm_index.size();

            left_rev = left;
            right_rev = right;
        } else {
            left = ind_start;
            right = ind_start + 1;

            left_rev = right;
            right_rev = left;
        }
    } else {
        left = num_begin;
        if (marker_value + 2 <= maxx)
            right = fm_index.C[fm_index.char2comp[marker_value + 2]];
        else
            right = fm_index.size();

        left_rev = left;
        right_rev = right;
    }
    return last;
}
