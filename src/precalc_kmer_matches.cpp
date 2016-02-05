#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>
#include <string>
#include <list>
#include <utility>
#include <boost/functional/hash.hpp>


using namespace sdsl;

//what to do with Ns?

void generate_all_kmers(std::vector<uint8_t> letters, std::vector<uint8_t>& substr, int k, int n, std::vector<std::vector<uint8_t>>& kmers);

void precalc_kmer_matches (csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa, int k,   
			   sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx, 
			   sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx_rev,
			   sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>& kmer_sites,
			   std::vector<int> mask_a, uint64_t maxx, sequence_set<std::vector<uint8_t>>& kmers_in_ref) 
{
  std::vector<uint8_t> letters; // add N/other symbols?
  std::vector<std::vector<uint8_t>> kmers;
  std::vector<uint8_t> substr;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> temp2;		  
  std::list<std::pair<uint64_t,uint64_t>> temp;
  bool first_del;

  for (uint8_t i=1;i<=4;i++) letters.push_back(i);
  generate_all_kmers(letters, substr,k, letters.size(),kmers);
  
  for (std::vector<std::vector<uint8_t>>::iterator it=kmers.begin();  it!=kmers.end(); ++it) {
    kmer_idx[(*it)]=temp;
    kmer_idx_rev[(*it)]=temp;
    kmer_sites[(*it)]=temp2;
    first_del=false;
    std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size()-1,0,csa.size()-1,(*it).begin(),(*it).end(),kmer_idx[(*it)],kmer_idx_rev[(*it)],kmer_sites[(*it)],mask_a,maxx,first_del);
    if (first_del) kmers_in_ref.insert(*it);
  }
}  




void generate_all_kmers(std::vector<uint8_t> letters, std::vector<uint8_t>& substr, int k, int n, std::vector<std::vector<uint8_t>>& kmers)
{
  if (k==1) {
    for (int j = 0; j < n; j++)
      substr.push_back(letters[j]);
      kmers.push_back(substr);
  }
  else {
    for (int jj = 0; jj < n; jj++)
      substr.push_back(letters[jj]);
      generate_all_kmers(letters, substr, k-1, n,kmers);
  }
}
