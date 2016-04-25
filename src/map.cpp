#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <fstream>
#include <iostream>
#include <cstdint>
#include <time.h> 
#include <vector>  
#include "bwt_search.h"
#include <seqread.hpp>
#include "precalc_gen.hpp"

using namespace std;  
using namespace sdsl;

void timestamp(); 

//argv[1] - file containing linear prg
//argv[2] - file where CSA is stored 
//argv[3] - file containing reads to be mapped (one read per line)
//argv[4] - file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which site
//argv[5] - file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which allele
//argv[6] - name of output file where coverages on each allele are printed
//argv[7] - name of output file where reads that have been processed are printed
//argv[8] - name of binary file where the prg in integer alphabet is written
//argv[9] - memory log file for CSA
//argv[10] - size of precalculated kmers
//argv[11] - kmer file

int main(int argc, char* argv[]) {

	std::vector<uint64_t> mask_s;
	std::vector<int> mask_a;
	std::vector<std::vector<int> > covgs;
	std::string q;

	SeqRead inputReads(argv[3]); 
	std::ofstream out(argv[6]); 
	std::ofstream out2(argv[7]);

	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>> kmer_idx,kmer_idx_rev;
	sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> kmer_sites;
	sequence_set<std::vector<uint8_t>> kmers_in_ref;

	//not using mask_s anymore?
	uint64_t maxx=parse_masks(mask_s,mask_a,argv[4],argv[5],covgs);

	timestamp();
	cout<<"CSA construction"<<endl;
	csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(argv[1],covgs,argv[8],argv[9],argv[2],true);
	timestamp();

	std::vector<std::vector<string> > site_reads(covgs.size(),std::vector<string>(1));
	int no_reads=0;
	uint64_t no_mapped=0;
	int inc=0;
	int no_occ=0;
	bool first_del=false;

	int k=atoi(argv[10]); //verify input
	get_precalc_kmers(csa,kmer_idx,kmer_idx_rev,kmer_sites,kmers_in_ref,mask_a,argv[11],maxx,k);

	//precalc_kmer_matches(csa,k,kmer_idx,kmer_idx_rev,kmer_sites,mask_a,maxx,kmers_in_ref,argv[11]);
	timestamp();

	std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
	std::list<std::pair<uint64_t,uint64_t>>::iterator it, it_rev;
	std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
	std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>::iterator it_s;
	std::vector<uint8_t>::iterator res_it;

	for (auto q: inputReads)
	{
		// If you declare p inside the while scope, it is destroyed/created automatically in every loop
		std::vector<uint8_t> p;

		//logging
		if (!(inc++%10)) { out2<<no_reads<<endl; }

		//add N's

		for (int i=0,seqlen=strlen(q->seq);i<seqlen;i++) {
			if (q->seq[i]=='A' or q->seq[i]=='a') p.push_back(1);
			if (q->seq[i]=='C' or q->seq[i]=='c') p.push_back(2);
			if (q->seq[i]=='G' or q->seq[i]=='g') p.push_back(3);
			if (q->seq[i]=='T' or q->seq[i]=='t') p.push_back(4);
		}

		std::vector<uint8_t> kmer(p.begin()+p.size()-k,p.end()); //is there a way to avoid making this copy?
		if (kmer_idx.find(kmer)!=kmer_idx.end() && kmer_idx_rev.find(kmer)!=kmer_idx_rev.end() && kmer_sites.find(kmer)!=kmer_sites.end()) {
		  sa_intervals=kmer_idx[kmer];
		  sa_intervals_rev=kmer_idx_rev[kmer];
		  sites=kmer_sites[kmer];	
		  
		  it=sa_intervals.begin();
		  it_rev=sa_intervals_rev.begin();

		  if (kmers_in_ref.find(kmer)!=kmers_in_ref.end()) first_del=false;
		  else first_del=true;

		  res_it=bidir_search_bwd(csa, (*it).first, (*it).second, (*it_rev).first, (*it_rev).second, p.begin(),p.begin()+p.size()-k, sa_intervals, sa_intervals_rev, sites, mask_a, maxx, first_del);

		  no_occ=0;
		  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
		    no_occ+=(*it).second-(*it).first;

		  sa_intervals.clear();
		  sa_intervals_rev.clear();
		  sites.clear();
		}
		else no_occ=0;
	   
		out<<no_occ<<" ";
		//clear p, sa_intervals etc
		p.clear();

		no_reads++;
	}
	/*
	for (int i=0;i<covgs.size();i++) {
		for (int j=0;j<covgs[i].size();j++)
			out<<covgs[i][j]<<" ";
		out<<endl;
	}

	for (int i=0;i<site_reads.size();i++) {
		std::string name=string(argv[9])+"Site"+std::to_string(i);
		std::ofstream fsite(name);
		for (int j=0;j<site_reads[i].size();j++)
			fsite<<site_reads[i][j]<<endl;
		fsite.close();
		}*/
	timestamp();

	out.close();
	out2.close();
	return(0);
}

void timestamp(){
	time_t ltime;
	ltime = time(NULL);
	printf("\n-----\n%s",asctime(localtime(&ltime)));
	fflush(stdout);
}
