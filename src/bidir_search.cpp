#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <cassert>
#include <tuple>
#include <cstdint>

#include "fm_index.hpp"
#include "ranks.hpp"


//csa is the compressed suffix array object 
//    for the text you're searching on, 
//    (based on a wavelet tree (WT) over the BWT)
//[left, right] is the SA interval of the occurences of the pattern 
//              you're trying to extend
//[left_rev, 
//  right_rev]  is the SA interval of the occurences of the pattern
//              you're trying to extend...
//              in the reverse text csa 
//              (don't need the actual reverse text csa, 
//              just these indices)
//char next_char for extending the current pattern


std::pair<uint64_t, uint64_t> bidir_search(const uint8_t next_char,
                                           const std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                                           const std::list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                                           const FM_Index &fm_index,
                                           const DNA_Rank &rank_all) {

    uint64_t left = sa_interval_it->first;
    uint64_t right = sa_interval_it->second;

    assert(left < right);
    assert(right <= fm_index.size());

    // next_char_interval_left (below) is the first occurrence/posn
    //          of char next_char in the far left column
    //          of the BW matrix
    const uint64_t next_char_interval_left = fm_index.C[fm_index.char2comp[next_char]];

    // [ NB Since the suffixes are alphabetically ordered,
    // the position at which next_char appears for the first time
    // in this first column is equal to the number of
    // times characters smaller than next_char appear in text  ]

    auto tmp = rank_all.find(next_char - 1)->second;
    if (left == 0)
        left = next_char_interval_left;
    else
        left = next_char_interval_left + tmp[left - 1];

    right = next_char_interval_left + tmp[right - 1];
    assert(right >= left);

    return std::make_pair(left, right);
}
