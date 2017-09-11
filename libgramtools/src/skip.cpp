#include <sdsl/suffix_arrays.hpp>

#include "prg.hpp"
#include "fm_index.hpp"
#include "skip.hpp"


bool skip(uint64_t &left, uint64_t &right, const uint64_t marker_value, const PRG_Info &prg_info) {

    const auto is_edge_marker = marker_value % 2 == 1;
    if (is_edge_marker) {
        bool last = process_variant_edge_marker(left, right, marker_value, prg_info);
        return last;
    }

    const uint64_t num_begin = prg_info.fm_index.C[prg_info.fm_index.char2comp[marker_value - 1]];
    const auto site_boundary_lower = prg_info.fm_index[num_begin];
    const auto site_boundary_upper = prg_info.fm_index[num_begin + 1];

    if (site_boundary_lower < site_boundary_upper) {
        left = num_begin;
        right = num_begin + 1;
    } else {
        left = num_begin + 1;
        right = num_begin + 2;
    }

    return false;
}

bool process_variant_edge_marker(uint64_t &left, uint64_t &right,
                                 const uint64_t marker_value,
                                 const PRG_Info &prg_info) {

    const uint64_t num_begin = prg_info.fm_index.C[prg_info.fm_index.char2comp[marker_value]];
    const auto site_boundary_lower = prg_info.fm_index[num_begin];
    const auto site_boundary_upper = prg_info.fm_index[num_begin + 1];

    uint64_t ind_start;
    if (site_boundary_lower > site_boundary_upper)
        ind_start = num_begin + 1;
    else
        ind_start = num_begin;

    bool last = false;
    if (right - left == 1) {
        const uint64_t site_end = std::max(site_boundary_lower, site_boundary_upper);
        if (prg_info.fm_index[left] == site_end + 1) {
            last = true;

            left = num_begin;
            if (marker_value + 2 <= prg_info.max_alphabet_num)
                right = prg_info.fm_index.C[prg_info.fm_index.char2comp[marker_value + 2]];
            else
                right = prg_info.fm_index.size();
        } else {
            left = ind_start;
            right = ind_start + 1;
        }
    } else {
        left = num_begin;
        if (marker_value + 2 <= prg_info.max_alphabet_num)
            right = prg_info.fm_index.C[prg_info.fm_index.char2comp[marker_value + 2]];
        else
            right = prg_info.fm_index.size();
    }
    return last;
}
