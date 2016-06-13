#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
//#include "bwt_search.h"
#include <tuple>
#include <cstdint>


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
//char c for extending the current pattern     
using namespace sdsl;

uint64_t bidir_search(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> &csa, 
		      uint64_t& left, uint64_t& right, 
		      uint64_t& left_rev, uint64_t& right_rev, 
		      uint8_t c)
{
  assert(left < right); 
  assert(right <= csa.size());
  assert((c>0) & (c<5));  //would be nice to replace 5 with a constant set at compile-time (so one day can do with amino); the n parameter in precalc_kmer_matches

  // c_begin (below) is the first occurrence/posn 
  //          of char c in the far left column 
  //          of the BW matrix 
  uint64_t c_begin = csa.C[csa.char2comp[c]];

  // [ NB Since the suffixes are alphabetically ordered, 
  // the position at which c appears for the first time 
  // in this first column is equal to the number of 
  // times characters smaller than c appear in text  ]

  //r_s_b is a tuple  - (rank_l,s,b)
  auto r_s_b =  csa.wavelet_tree.lex_count(left, right, c);

  // r stands for rank
  // s - characters Smaller than c, 
  // b - characters Bigger than c


  //rank_l is the number of times c occurs in BWT[0,left]
  uint64_t rank_l = std::get<0>(r_s_b);   
  // s is total number of times 
  //      chars smaller than c occur in BWT[left,right]
  uint64_t s = std::get<1>(r_s_b);
  // b is total number of times 
  //      chars bigger than c occur in BWT[left,right]  
  uint64_t b = std::get<2>(r_s_b);

  //rank_r is the number of times c occurs in BWT[0,right]
  uint64_t rank_r = right - left - s - b + rank_l - 1;

  //get new interval, after appending c
  left  = c_begin + rank_l;
  right = c_begin + rank_r + 1;
  assert(right>=left);

  //now same in reverse csa
  left_rev  = left_rev + s;
  right_rev = right_rev - b + 1;
  //  assert(right_rev-left_rev == right-left);

  return right-left;
}
