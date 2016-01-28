#ifndef __BWT_SEARCH_H_INCLUDED__
#define __BWT_SEARCH_H_INCLUDED__

#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <boost/functional/hash.hpp> 

using namespace sdsl;

template < typename SEQUENCE > struct seq_hash
			       {
				 std::size_t operator() ( const SEQUENCE& seq ) const
				 {
				   std::size_t hash = 0 ;
				   boost::hash_range( hash, seq.begin(), seq.end() ) ;
				   return hash ;
				 }
};

template < typename SEQUENCE, typename T >

using sequence_map = std::unordered_map< SEQUENCE, T, seq_hash<SEQUENCE> > ;

csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_constr(std::string fname, std::vector<std::vector<int>>& covgs, char* int_al_fname, char* memory_log_fname, char* csa_file);

void parse_masks(std::vector<uint64_t>& mask_s, std::vector<int>& mask_a, std::string sites_fname, std::string alleles_fname, std::vector<std::vector<int>>& covgs);

uint64_t bidir_search(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa, 
		      uint64_t& left, uint64_t& right, 
		      uint64_t& left_rev, uint64_t& right_rev, 
		      uint8_t c);


std::pair<uint32_t, std::vector<int>> get_location(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa,
		      uint64_t num_idx,
		      uint32_t num, bool last,
		      std::vector<int>& allele,
		      std::vector<int> mask_a);

bool skip(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa,
                      uint64_t& left, uint64_t& right,
                      uint64_t& left_rev, uint64_t& right_rev,
	              uint32_t num);

std::vector<uint8_t>::iterator bidir_search_bwd(csa_wt<wt_int<bit_vector,rank_support_v5<>>> csa,
						uint64_t left, uint64_t right,
						uint64_t left_rev, uint64_t right_rev,
						std::vector<uint8_t>::iterator pat_begin, std::vector<uint8_t>::iterator pat_end,
						std::list<std::pair<uint64_t,uint64_t>>& sa_intervals,
						std::list<std::pair<uint64_t,uint64_t>>& sa_intervals_rev,
						std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>& sites,
						std::vector<int> mask_a, uint64_t maxx);

#endif
