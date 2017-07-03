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
#include <fstream>

#include "map.hpp"
#include "bwt_search.h"
#include "definitions.hpp"


using namespace sdsl;

//what to do with Ns?

//void generate_all_kmers(std::vector<uint8_t> letters, std::vector<uint8_t>& substr, int k, int n, std::vector<std::vector<uint8_t>>& kmers);
//void get_kmers(char *kmerfile, std::vector<std::vector<uint8_t>>& kmers);

void precalc_kmer_matches (const csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> &csa, int k,
		sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx, 
		sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx_rev,
		sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>& kmer_sites,
		std::vector<int> &mask_a, uint64_t maxx, sequence_set<std::vector<uint8_t>>& kmers_in_ref, std::vector<std::vector<uint8_t>> &kmers,
        const VariantMarkers &variants, int thread_id) {

	std::cout << std::setprecision(2) << std::fixed;
	std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> temp2;		  
	std::list<std::pair<uint64_t,uint64_t>> temp;
	bool first_del;

	for (unsigned long i=0; i < kmers.size(); i++) {
		auto &kmer = kmers[i];

		if (thread_id == LOG_THREAD_ID) {
            std::cout << thread_id << ": kmers processed: "
                      << i << "/" << kmers.size() << std::endl;
        }

		kmer_idx[kmer]=temp;
		kmer_idx_rev[kmer]=temp;
		kmer_sites[kmer]=temp2;
		first_del=false;
		bool precalc_done=false;

		bidir_search_bwd(csa,0,
                         csa.size(),0,
                         csa.size(),
                         (kmer).begin(),(kmer).end(),
                         kmer_idx[kmer],kmer_idx_rev[kmer],
                         kmer_sites[kmer], mask_a,maxx,first_del,
                         precalc_done, variants);
		if  ((kmer_idx[kmer]).empty())
		  {
		    kmer_idx.erase(kmer);
		  }
		if  ((kmer_idx_rev[kmer]).empty())
		  {
		    kmer_idx_rev.erase(kmer);
		  }
		
		if (!first_del) kmers_in_ref.insert(kmer);
		
	}
}
