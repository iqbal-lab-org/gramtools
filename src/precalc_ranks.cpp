#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include <cstdint>

using namespace sdsl;
using namespace std;


void precalc_ranks(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216>& csa,
  unordered_map<uint8_t,vector<uint64_t>>& rank_all)
{
  uint64_t bwt_size=csa.size();
  vector<uint64_t> rank(4,0);
  uint8_t symbols[]={1,2,3,4};

  for (auto c : symbols)
    rank_all[c-1]=vector<uint64_t> (bwt_size,0);

  for (auto i=0; i<bwt_size;i++) {
    auto curr_c=csa.bwt[i];
    if ((curr_c>0) && (curr_c<5)) {
      rank[curr_c-1]+=1;
      for (auto c:symbols) 
        rank_all[c-1][i]=rank[c-1];
    }
    else {
      for (auto c:symbols) 
        rank_all[c-1][i]=rank_all[c-1][i-1];
    }
  }
}
