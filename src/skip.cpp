#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>
#include <algorithm>

using namespace sdsl;

uint64_t skip(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa,
                      uint64_t& left, uint64_t& right,
                      uint64_t& left_rev, uint64_t& right_rev,
                      uint32_t num)
{
  uint64_t site_end,site_start;
  bool first;

  assert(left < right);
  assert(right <= csa.size());
  assert(num > 4);

  if (num%2==1) {
    uint64_t num_begin = csa.C[csa.char2comp[num]];

    site_end=max(csa[num_begin],csa[num_begin+1]); 
    site_start=min(csa[num_begin],csa[num_begin+1]);

    if (right-left==1) {
      if (csa[i]==csa[site_end]+1) {
	last=T

	left=num_begin;
	right=csa.C[csa.char2comp[num+2]];

	left_rev=left;
	right_rev=right;
      }
      else {
	//site=num;
	//allele.push_back(1);

	left=site_start;
	right=site_start+1;

	left_rev=right;
	right_rev=left;
	//	first=T;
      }
    }
    else {
      //site=num;
      //allele.push_back(1);

      left=num_begin;
      right=csa.C[csa.char2comp[num+2]];      

      left_rev=left;
      right_rev=right;
      //first=T;
    }
  }
  else {
    uint64_t num_begin = csa.C[csa.char2comp[num-1]];
    
    site_start=min(csa[num_begin],csa[num_begin+1]);
    site_end=max(csa[num_begin],csa[num_begin+1]);
    /*
    if (right-left==1) {
      site=num-1;
      allele.push_back(mask_a[csa[left]]);
    }
    else {
      site=num-1;
      for (i=left;i<right;i++)
	allele.push_back(mask_a[csa[i]]);
    }
    */
    left=site_start;
    right=site_start+1;

    left_rev=site_end;
    right_rev=site_end+1;
  }

  assert(right>=left);
  assert(right_rev-left_rev == right-left);

  return (last);
}
