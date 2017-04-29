#ifndef GRAMTOOLS_KMERS_HPP
#define GRAMTOOLS_KMERS_HPP

typedef sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>> KmerIdx;
typedef sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> KmerSites;
typedef sequence_set<std::vector<uint8_t>> KmersRef;


inline bool fexists (const std::string& name);

//trim from start
static inline std::string &ltrim(std::string &s);

// trim from end
static inline std::string &rtrim(std::string &s);

// trim from both ends
static inline std::string &trim(std::string &s);

std::vector<std::string> split(std::string cad,std::string delim);

using namespace std;  

#define THREADS 25

struct thread_data{
	csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> *csa;
	int k;
	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>* kmer_idx;
	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>* kmer_idx_rev;
	sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>* kmer_sites;
	std::vector<int> *mask_a ;
	uint64_t maxx ;
	sequence_set<std::vector<uint8_t>>*kmers_in_ref ;
	std::vector<std::vector<uint8_t>>*kmers;
};

void * worker (void *st);


void gen_precalc_kmers(
		csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> &csa,
//		sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx,
//		sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx_rev,
//		sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>& kmer_sites,
//		sequence_set<std::vector<uint8_t>> &kmers_in_ref,
//		std::vector<std::vector<uint8_t>> &kmers,
		std::vector<int> &mask_a,
		std::string kmer_fname,
		uint64_t maxx,
		int k
		);

void read_precalc_kmers(std::string fil, sequence_map<std::vector<uint8_t>, 
			std::list<std::pair<uint64_t,uint64_t>>> &kmer_idx, 
			sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>> &kmer_idx_rev, 
			sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, 
			std::vector<int>>>>> &kmer_sites,
			sequence_set<std::vector<uint8_t>>&kmers_in_ref
	);

void get_precalc_kmers(
		csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> &csa,
		sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx,
		sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx_rev,
		sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>& kmer_sites,
		sequence_set<std::vector<uint8_t>> &kmers_in_ref,
		std::vector<int> &mask_a,
		std::string kmer_fname,
		uint64_t maxx,
		int k);

#endif //GRAMTOOLS_KMERS_HPP
