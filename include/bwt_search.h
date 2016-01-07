#ifndef __BWT_SEARCH_H_INCLUDED__
#define __BWT_SEARCH_H_INCLUDED__

#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"

using namespace sdsl;

uint64_t bidir_search(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa, 
		      uint64_t& left, uint64_t& right, 
		      uint64_t& left_rev, uint64_t& right_rev, 
		      uint8_t c);


std::pair<uint32_t, std::vector<int>> get_location(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa,
		      uint64_t num_idx,
		      uint32_t num, bool last,
		      std::vector<int>& allele,
		      std::vector<int> mask_a);

bool skip(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa,
                      uint64_t& left, uint64_t& right,
                      uint64_t& left_rev, uint64_t& right_rev,
	              uint32_t num);

#endif
