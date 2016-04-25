#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <fstream>
#include <iostream>
#include <cstdint>
#include <time.h> 
#include <vector>  
#include "bwt_search.h"
#include <seqread.hpp>
#include <pthread.h>

using namespace std;  
using namespace sdsl;

#define THREADS 25

void timestamp(); 

//argv[1] - file containing linear prg
//argv[2] - file where CSA is stored 
//argv[4] - file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which site
//argv[5] - file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which allele
//argv[8] - name of binary file where the prg in integer alphabet is written
//argv[9] - memory log file for CSA
//argv[10] - size of precalculated kmers
//argv[11] - kmer file
//

struct thread_data{
	csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> *csa;
	int k;
	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>* kmer_idx;
	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>* kmer_idx_rev;
	sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>* kmer_sites;
	std::vector<int> *mask_a ;
	uint64_t maxx ;
	sequence_set<std::vector<uint8_t>>*kmers_in_ref ;
	std::vector<std::vector<uint8_t>>*kmers;
};

void * worker (void *st)
{
	thread_data* th=(thread_data *) st;
	precalc_kmer_matches(*(th->csa),th->k,*(th->kmer_idx),*(th->kmer_idx_rev),*(th->kmer_sites),*(th->mask_a),th->maxx,*(th->kmers_in_ref),*(th->kmers));
}

int main(int argc, char* argv[]) {
	int i=0;

	std::vector<uint64_t> mask_s;
	std::vector<int> mask_a;
	std::vector<std::vector<int> > covgs;
	pthread_t threads[THREADS];
	struct thread_data td[THREADS];


	std::vector<std::vector<uint8_t>> kmers[THREADS];

	ifstream kfile;
	string line;
	kfile.open(argv[11]);
	while (std::getline(kfile,line))
	{
		std::vector<uint8_t> kmer;
		for (auto c: line)
			switch (c)
			{
				case 'A': case 'a': kmer.push_back(1);break;
				case 'C': case 'c': kmer.push_back(2);break;
				case 'G': case 'g': kmer.push_back(3);break;
				case 'T': case 't': kmer.push_back(4);break;
			}
		kmers[i++].push_back(kmer);
		i%=THREADS;
	}

	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>> kmer_idx[THREADS],kmer_idx_rev[THREADS];
	sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> kmer_sites[THREADS];
	sequence_set<std::vector<uint8_t>> kmers_in_ref[THREADS];

	//not using mask_s anymore?
	uint64_t maxx=parse_masks(mask_s,mask_a,argv[4],argv[5],covgs);

	cout<<"CSA construction"<<endl;
	csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(argv[1],covgs,argv[8],argv[9],argv[2],true);

	int k=atoi(argv[10]); //verify input
	for (int i=0;i<THREADS;i++)
	{
		td[i].csa=&csa;
		td[i].k=k;
		td[i].kmer_idx=&kmer_idx[i];
		td[i].kmer_idx_rev=&kmer_idx_rev[i];
		td[i].kmer_sites=&kmer_sites[i];
		td[i].mask_a=&mask_a;
		td[i].maxx=maxx;
		td[i].kmers_in_ref=&kmers_in_ref[i];
		td[i].kmers=&kmers[i];
		pthread_create(&threads[i], NULL, worker, &td[i] );
//		worker(&td[i]);
//		thread_data* th=td+i;
//		precalc_kmer_matches(*(th->csa),th->k,*(th->kmer_idx),*(th->kmer_idx_rev),*(th->kmer_sites),*(th->mask_a),th->maxx,*(th->kmers_in_ref),*(th->kmers));
		std::cout << i << std::endl;
//		pthread_join(threads[i],&status);
	}

	for (int i=0;i<THREADS;i++){
		void * status;
		pthread_join(threads[i],&status);
		for (auto k: kmers[i])
		{
			for (auto n: k) std::cout<<(int) n << ' ';
			for (auto o: kmer_idx[i][k]) std::cout << o.first << ' ' << o.second << ' ';
			std::cout << '|';
			for (auto o: kmer_idx_rev[i][k]) std::cout << o.first << ' ' << o.second << ' ';
			std::cout << '|';
			for (auto o: kmer_sites[i][k]) 
			{
				for (auto v: o)
				{ 
					std::cout << v.first << ' ';
					for (auto n: v.second) std::cout<<(int) n << ' ';
					std::cout << '@';
				}
				std::cout << '|';
			}

			std::cout << std::endl;
		}
	}


	return(0);
}

void timestamp(){
	time_t ltime;
	ltime = time(NULL);
	printf("\n-----\n%s",asctime(localtime(&ltime)));
	fflush(stdout);
}
