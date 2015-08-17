#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>

using namespace sdsl;

std::vector<uint8_t>::iterator bidir_search_bwd(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa,
                      uint64_t& left, uint64_t& right,
		      uint64_t& left_rev, uint64_t& right_rev,
                      std::vector<uint8_t>::iterator pat_begin, std::vector<uint8_t>::iterator pat_end)
{
  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev; 
  std::vector<uint8_t>::iterator pat_it=pat_end;
  std::list<int>::iterator it2, it2_rev;
  uint8_t c=*pat_it;
  bool last, first_del=F;

  sa_intervals.push_back(std::make_pair(left,right));
  sa_intervals_rev.push_back(std::make_pair(left_rev,right_rev);

  while (pat_it>pat_begin && !sa_intervals.empty()) {
    --pat_it;
    
    assert(sa_intervals.size()==sa_intervals_rev.size());
    for(it=sa_intervals.begin(), it_rev=sa_intervals_rev.begin(); it!=sa_intervals.end() && it_rev!=sa_intervals_rev.end(); ++it, ++it_rev) {
      auto res= csa.wavelet_tree.range_search_2d(*it.first, *it.second-1, 5, max).second;
      //res is going to contain sometimes multiple 5's for ex, non adjacent
      for (auto z : res) {
	auto i=z.first;
	auto num=z.second;

	left_new=*it.first;
	right_new=*it.second;

	left_rev_new=*it_rev.first;
	right_rev_new=*it_rev.second;

	last=skip(csa,left_new,right_new,left_rev_new,right_rev_new,num);
	//call location here
	//need original [l,r] to for the next loop iterations
	// how to alternate between forward and backward?
	if (it==sa_intervals.begin() && first_del==F) {
	  sa_intervals.push_back(std::make_pair(left_new,right_new));
	  sa_intervals_rev.push_back(std::make_pair(left_rev_new,right_rev_new));
	}
	else {
	  *it=std::make_pair(left_new,right_new);
	  *it_rev=std::make_pair(left_rev_new,right_rev_new);
	}
      }
    }

    assert(sa_intervals.size()==sa_intervals_rev.size());
    it=sa_intervals.begin();
    it_rev=sa_intervals_rev.begin();			
		  
    while (it!=sa_intervals.end() && it2_rev!=sa_intervals_rev.end()) {	
      //calculate sum to return
      if (bidir_search(csa,*it.first,*it.second,*it_rev.first,*it_rev.second,c)>0) {
	++it2;
	++it2_rev;
      }
      else {
	if (it==sa_intervals.begin()) first_del=T;
	it=sa_intervals.erase(it);
	it_rev=sa_intervals_rev.erase(it);
      }
    }
  }
  
  if (pat_it!=pat_begin) return(pat_it); // where it got stuck
  else {
    if !sa_intervals.empty() return(pat_end);
    else return(pat_begin); //where it got stuck
  }
			     }			     			          
