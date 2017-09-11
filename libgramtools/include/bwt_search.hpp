#include <cstdlib>
#include <iostream>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "kmers.hpp"
#include "ranks.hpp"


#ifndef GRAMTOOLS_BWT_SEARCH_H
#define GRAMTOOLS_BWT_SEARCH_H

void precalc_ranks(const FM_Index &fm_index, const DNA_Rank &rank_all);

SA_Interval bidir_search(const uint8_t next_char, const SA_Interval &sa_interval, const PRG_Info &prg_info);

bool skip(uint64_t &left, uint64_t &right, const uint64_t num, const PRG_Info &prg_info);


#endif //GRAMTOOLS_BWT_SEARCH_H
