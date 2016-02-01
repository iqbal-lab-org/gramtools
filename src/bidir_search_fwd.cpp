#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>

using namespace sdsl;

std::vector<uint8_t>::iterator bidir_search_bwd(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa_rev,
                      uint64_t left, uint64_t right,
		      uint64_t left_rev, uint64_t right_rev,
		      std::vector<uint8_t>::iterator pat_begin, std::vector<uint8_t>::iterator pat_end,
		      std::list<std::pair<uint64_t,uint64_t>>& sa_intervals, 
		      std::list<std::pair<uint64_t,uint64_t>>& sa_intervals_rev,
		      std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>& sites,
		      std::vector<int> mask_a, uint64_t maxx, bool& first_del)

//need to swap * with *_rev everywhere
{
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>::iterator it_s;
  std::vector<uint8_t>::iterator pat_it=pat_begin;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it, it_rev;
  uint8_t c=*pat_it;
  bool last,ignore;
  uint64_t left_new, right_new, left_rev_new, right_rev_new;
  std::vector<std::pair<uint32_t, std::vector<int>>> empty_pair_vector;
  std::vector<int> allele_empty;
  std::vector<std::pair<uint64_t,uint64_t>> res;  

  assert(left<right);
  assert(right<=csa_rev.size());

  sa_intervals.push_back(std::make_pair(left,right));
  sa_intervals_rev.push_back(std::make_pair(left_rev,right_rev));
  sites.push_back(empty_pair_vector);

  while (pat_it<pat_end && !sa_intervals.empty()) {

    assert(sa_intervals.size()==sa_intervals_rev.size());
    assert(sa_intervals.size()==sites.size());

    it=sa_intervals.begin();
    it_rev=sa_intervals_rev.begin();
    it_s=sites.begin();

    for(;it!=sa_intervals.end() && it_rev!=sa_intervals_rev.end() && it_s!=sites.end(); ++it, ++it_rev, ++it_s) {
      res= csa_rev.wavelet_tree.range_search_2d((*it).first, (*it).second-1, 5, maxx).second;
      //might want to sort res based on pair.second - from some examples it looks like sdsl already does that so res is already sorted 
      uint32_t prev_num=0;
      for (auto z : res) {
	uint64_t i=z.first;
	uint32_t num=z.second;

	if (num==prev_num && num%2==0) ignore=true;
	else ignore=false;

	left_new=(*it).first;
	right_new=(*it).second;

	//need original [l,r] to for the next loop iterations

	left_rev_new=(*it_rev).first;
	right_rev_new=(*it_rev).second;

	last=skip(csa_rev,left_new,right_new,left_rev_new,right_rev_new,num);
	
	// how to alternate between forward and backward?
	if (it==sa_intervals.begin() && first_del==false && !ignore) {
	  sa_intervals.push_back(std::make_pair(left_new,right_new));
	  sa_intervals_rev.push_back(std::make_pair(left_rev_new,right_rev_new));
	  sites.push_back(std::vector<std::pair<uint32_t, std::vector<int>>>(1,get_location(csa_rev,i,num,last,allele_empty,mask_a)));
	}
	    //there will be entries with pair.second empty (corresp to allele) coming from crossing the last marker
	    //can delete them here or in top a fcn when calculating coverages
	else {
	  if (ignore) sites.back().back()=get_location(csa_rev,i,num,last,sites.back().back().second,mask_a);
	  else {
	    *it=std::make_pair(left_new,right_new);
	    *it_rev=std::make_pair(left_rev_new,right_rev_new);
	    (*it_s).push_back(get_location(csa_rev,i,num,last,allele_empty,mask_a));
	  }
	}
	prev_num=num;  
      }
     }
    

    assert(sa_intervals.size()==sa_intervals_rev.size());
    assert(sa_intervals.size()==sites.size());

    it=sa_intervals.begin();
    it_rev=sa_intervals_rev.begin();	
    it_s=sites.begin();
		  
    while (it!=sa_intervals.end() && it_rev!=sa_intervals_rev.end() && it_s!=sites.end()) {	
      //calculate sum to return- can do this in top fcns
      if (bidir_search(csa_rev,(*it).first,(*it).second,(*it_rev).first,(*it_rev).second,c)>0) {
	++it;
	++it_rev;
	++it_s;
      }
      else {
	if (it==sa_intervals.begin()) first_del=true;
	//might need to see first_del from top fcns to check if there are matches in the reference
	it=sa_intervals.erase(it);
	it_rev=sa_intervals_rev.erase(it_rev);
	it_s=sites.erase(it_s);
      }
    }
    ++pat_it;
  }  

  if (pat_it!=pat_begin) return(pat_it); // where it got stuck
  else {
    if (!sa_intervals.empty()) return(pat_end);
    else return(pat_begin); //where it got stuck
  }
  
}			     			     			          
