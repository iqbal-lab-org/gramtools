#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>

using namespace sdsl;

std::pair<uint32_t, std::vector<int>> get_location(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa,
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

  //problems:
  //avoid mask_s
  //when have matches to the right of both odd numbers
