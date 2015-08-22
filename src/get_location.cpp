#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
//#include "bwt_search.h"
#include <tuple>
#include <cstdint>

using namespace sdsl;

std::pair<uint32_t, std::vector<int>> get_location(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa,
		      uint64_t num_idx,
		      uint32_t num, bool last,
		      std::vector<int>& allele,
		      std::vector<int> mask_a)
{
  uint32_t site;

  if (num%2==1) {
 	site=num;
        if (!last) allele.push_back(1);
      }
  else {
      site=num-1;
      allele.push_back(mask_a[csa[num_idx]]);
    }
  return (std::make_pair(site,allele));    
}  

//removed left and right args because they only form a proper interval is when last is true; but then they can't help with location because allele is unknown and we're just addin site
//when allele is added, the non-last markers are encountered from the left and at that point the interval is either of length 1 (if we're comeng from the left of the last num), or we need to split it in 1-length subintervals to add each allele; might be worth looking whether avoiding this splitting can make things more efficient

  //problems:
  //avoid mask_s
  //when have matches to the right of both odd numbers
  //when result is an interval, each element needs to have location saved - in current framework, result is never going to be an interval because res splits it
