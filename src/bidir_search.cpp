#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>


//csa is the compressed suffix array object for the text you're searching on, 
//    based on a wavelet tree (WT) over the BWT
//[left, right]         is the SA interval of the occurences of the pattern you're trying to extend
//[left_rev, right_rev] is the SA interval of the occurences of the pattern you're trying to extend 
//                      in the reverse text csa (don't need the actual reverse text csa, just these indices) 
//char c for extending the current pattern     

uint32_t bidir_search(csa_wt<wt_int<rrr_vector<63>>> csa, 
		      uint32_t& l, uint32_t& r, 
		      uint32_t& l_rev, uint32_t& r_rev, 
		      uint32_t c)
{
  assert(l < r); assert(r <= csa.size());


  // the beginning position of char c in the First column <<< ZAM says - what?
  uint32_t c_begin = csa.C[csa.char2comp[c]];//total number of occurences of chars smaller than c in text;  << ????????????????



  //r_s_-b is a tuple??<<<<<<
  auto r_s_b =  csa.wavelet_tree.lex_count(l, r, c);// tuple (r,s,b)
  uint32_t rank_l = std::get<0>(r_s_b);// r is number of times c occurs in BWT[l,r]
  uint32_t s = std::get<1>(r_s_b); // s is total number of times chars smaller than c occur in BWT[l,r]
  uint32_t b = std::get<2>(r_s_b);// b is total number of times chars bigger than c occur in BWT[l,r] 
  uint32_t rank_r = r - l - s - b + rank_l;
  l = c_begin + rank_l;
  r = c_begin + rank_r+1;
  assert(r >=l);
  l_rev = l_rev + s;
  r_rev = r_rev - b+1;
  assert(r_rev-l_rev == r-l);
  return r-l;
}
