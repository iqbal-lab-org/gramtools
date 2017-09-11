#include <sdsl/suffix_arrays.hpp>

#include "kmers.hpp"


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


SA_Interval bidir_search(const uint8_t next_char, const SA_Interval &sa_interval, const PRG_Info &prg_info) {

    uint64_t left = sa_interval.first;
    uint64_t right = sa_interval.second;

    assert(left < right);
    assert(right <= prg_info.fm_index.size());

    // next_char_interval_left (below) is the first occurrence/posn
    //          of char next_char in the far left column
    //          of the BW matrix
    const uint64_t next_char_interval_left = prg_info.fm_index.C[prg_info.fm_index.char2comp[next_char]];

    // [ NB Since the suffixes are alphabetically ordered,
    // the position at which next_char appears for the first time
    // in this first column is equal to the number of
    // times characters smaller than next_char appear in text  ]

    auto tmp = prg_info.dna_rank.find(next_char - 1)->second;
    if (left == 0)
        left = next_char_interval_left;
    else
        left = next_char_interval_left + tmp[left - 1];

    right = next_char_interval_left + tmp[right - 1];
    assert(right >= left);

    return std::make_pair(left, right);
}
