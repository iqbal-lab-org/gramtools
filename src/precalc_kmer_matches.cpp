#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>
#include <string>

using namespace sdsl;

//what to do with Ns?

void generate_all_kmers(std::vector<std::string> letters, std::string substr, int k, int n, std::vector<std::string>& kmers);

void precalc_kmer_matches (csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa, int k,   
			   std::unordered_map<std::string, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx, 
			   std::unordered_map<std::string, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx_rev,
			   std::unordered_map<std::string, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>& kmer_sites) 
{
  std::vector<std::string> letters={"1","2","3","4"}; // add N/other symbols?
  std::vector<std::string> kmers;

  generate_all_kmers(letters, "",k, letters.size(),kmers);
  
  for (std::vector<std::string>::iterator it=kmers.begin();  it!=kmers.end(); ++it) {
    kmer_idx[*it]=std::list<std::pair<uint64_t,uint64_t>> ();
    kmer_idx_rev[*it]=std::list<std::pair<uint64_t,uint64_t>> ();
    kmer_sites[*it]=std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> ();
    auto res_it=bidir_search_bwd(csa,0,csa.size()-1,0,csa.size()-1,(*it).begin(),(*it).end,kmer_idx[*it],kmer_idx_rev[*it],kmer_sites[*it]);
  }
}  




void generate_all_kmers(std::vector<std::string> letters, std::string substr, int k, int n, std::vector<std::string>& kmers)
{
  if (k==1) {
    for (int j = 0; j < n; j++)
      kmers.push_back(substr+letters[j]);
  }
  else {
    for (int jj = 0; jj < n; jj++)
      generate_all_kmers(letters, substr+letters[jj],k-1, n,kmers);
  }
}
