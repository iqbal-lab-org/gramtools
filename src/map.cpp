#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <time.h> 
#include <vector>  
#include "bwt_search.h"
#include <seqread.hpp>
#include "precalc_gen.hpp"
#include <getopt.h>

using namespace std;  
using namespace sdsl;

void timestamp(); 


char *HELP="\ngramtools usage:\n\
\n	  --prg   -p   input file containing linear prg                                                                                                     \
\n	  --csa   -c   output file where CSA is stored                                                                                                       \
\n	  --input -i   input FASTA/FASTQ file to be mapped                                                                         \
\n	  --ps    -s   input file containing mask over the                                                                        \      \n                        linear prg that indicates at each                                                                       \      \n                        position whether you are inside a                                                                       \      \n                        site and if so, which site   \
\n	  --pa    -a   input file containing mask over the                                                                        \      \n                        linear prg that indicates at each                                                                       \      \n                        position whether you are inside a                                                                       \      \n                        site and if so, which allele   \
\n	  --co    -v   name of output file where coverages on each allele are printed                                                                 \
\n	  --ro    -r   name of output file where reads that have been processed are printed                                                           \
\n	  --po    -b   output filename of binary file containing the prg in integer alphabet                                                                \
\n	  --log   -l   Output memory log file for CSA                                                                                                        \
\n	  --ksize -k   size of precalculated kmers                                                                                                    \
\n	  --kfile -f   input  file listing all kmers in PRG                                                                                                                     \
\n";


//argv[1] -  file containing linear prg
//argv[2] -  file where CSA is stored 
//argv[3] -  file containing reads to be mapped (one read per line)
//argv[4] -  file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which site
//argv[5] -  file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which allele
//argv[6] -  name of output file where coverages on each allele are printed
//argv[7] -  name of output file where reads that have been processed are printed
//argv[8] -  name of binary file where the prg in integer alphabet is written
//argv[9] -  memory log file for CSA
//argv[10] - size of precalculated kmers
//argv[11] - kmer file

int main(int argc, char* argv[]) {

	int c;

	std::string _prg="",_csa="",_input="",_sitemask="",_allmask="",_covoutput="",_readoutput="",_binoutput="",_log="",_ksize="",_kfile="";
	std::vector<std::string *>pars={&_prg,&_csa,&_input,&_sitemask,&_allmask,&_covoutput,&_readoutput,&_binoutput,&_log,&_ksize,&_kfile};

	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			//	{"verbose", no_argument,       &verbose_flag, 1},
			//	{"brief",   no_argument,       &verbose_flag, 0},
			/* These options don’t set a flag.
			 *              We distinguish them by their indices. */
			{"prg",     required_argument, 0, 'p'},
			{"csa",     required_argument, 0, 'c'},
			{"input",   required_argument, 0, 'i'},
			{"ps"  ,    required_argument, 0, 's'},
			{"pa"  ,    required_argument, 0, 'a'},
			{"co"  ,    required_argument, 0, 'v'},
			{"ro"  ,    required_argument, 0, 'r'},
			{"po"  ,    required_argument, 0, 'b'},
			{"log"  ,    required_argument, 0, 'l'},
			{"ksize"  ,    required_argument, 0, 'k'},
			{"kfile"  ,    required_argument, 0, 'f'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		c = getopt_long (argc, argv, "p:c:i:s:a:b:v:r:l:k:f:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
	//		case 0:
	//			/* If this option set a flag, do nothing else now. */
	//			if (long_options[option_index].flag != 0)
	//				break;
	//			printf ("option %s", long_options[option_index].name);
	//			if (optarg)
	//				printf (" with arg %s", optarg);
	//			printf ("\n");
	//			break;

			case 'p': case 'c': case 'i': case 's': case 'a': case 'b': case 'v': case 'r': case 'l': case 'k': case 'f':
				uint16_t i;
				for (i=0;i<pars.size();i++)
				{
					if (long_options[i].val==c)
						*pars[i]=optarg;
				}
				break;

			case '?':
				/* getopt_long already printed an error message. */
				std::cout << "Error parsing arguments\n" << HELP;
				break;

			default:
				std::cout << "Error parsing arguments\n" << HELP;
				abort ();
				break;
		}
	}

	for (auto i: pars)
		if (*i=="")
		{
			std::cout << "You must specify all parameters" << HELP;
			exit(-1);
		}

	std::vector<uint64_t> mask_s;
	std::vector<int> mask_a;

	std::vector<std::vector<float> > covgs;
	std::string q;

	SeqRead inputReads(_input.c_str()); 
	std::ofstream out(_covoutput); 
	std::ofstream out2(_readoutput);

	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>> kmer_idx,kmer_idx_rev;
	sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> kmer_sites;
	sequence_set<std::vector<uint8_t>> kmers_in_ref;

	//not using mask_s anymore?

	timestamp();
	cout<<"Start CSA construction"<<endl;
	csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=csa_constr(_prg,_binoutput,_log,_csa,true, true);
	timestamp();
	cout<<"End CSA construction"<<endl;

	uint64_t maxx=parse_masks(mask_s,mask_a,_sitemask,_allmask,covgs);

	std::vector<std::vector<string> > site_reads(covgs.size(),std::vector<string>(1));
	int no_reads=0;
	uint64_t no_mapped=0;
	int inc=0;
	int no_occ=0;
	bool invalid, first_del=false;

	int k=atoi(_ksize.c_str()); //verify input
	get_precalc_kmers(csa,kmer_idx,kmer_idx_rev,kmer_sites,kmers_in_ref,mask_a,_kfile,maxx,k);

	//precalc_kmer_matches(csa,k,kmer_idx,kmer_idx_rev,kmer_sites,mask_a,maxx,kmers_in_ref,_kfile);
	cout << "About to start mapping:"<<endl;
	timestamp();

	std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
	std::list<std::pair<uint64_t,uint64_t>>::iterator it, it_rev;
	std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
	std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>::iterator it_ss;
	std::unordered_set<int> repeats;
	std::vector<uint8_t>::iterator res_it;
	int in_sites;

	std::vector<uint8_t> p;//p is the read in integer alphabet
	p.reserve(200);	
	for (auto q: inputReads)
	{
	  
	  //logging
	  if (!(inc++%10)) { out2<<no_reads<<endl; }
	  
	  //add N's
	  int flag=0;

	  cout<<q->seq<<endl;
	  int seqlen=strlen(q->seq);
	  for (int i=0;i<seqlen;i++) {
	
	    if (q->seq[i]=='A' or q->seq[i]=='a') p.push_back(1);
	    else if (q->seq[i]=='C' or q->seq[i]=='c') p.push_back(2);
	    else if (q->seq[i]=='G' or q->seq[i]=='g') p.push_back(3);
	    else if (q->seq[i]=='T' or q->seq[i]=='t') p.push_back(4);
	    else {flag=1;}; 
	  }

	  if (flag==1)
	    {
	      continue;
	    }
	  std::vector<uint8_t> kmer(p.begin()+p.size()-k,p.end()); //is there a way to avoid making this copy?
	  if (kmer_idx.find(kmer)!=kmer_idx.end() && kmer_idx_rev.find(kmer)!=kmer_idx_rev.end() && kmer_sites.find(kmer)!=kmer_sites.end()) {
	          sa_intervals=kmer_idx[kmer];
		  sa_intervals_rev=kmer_idx_rev[kmer];
		  sites=kmer_sites[kmer];	
		  
		  it=sa_intervals.begin();
		  it_rev=sa_intervals_rev.begin();

		  //kmers in ref means kmers that do not cross any numbers 
		  //These are either in non-variable region, or are entirely within alleles
		  if (kmers_in_ref.find(kmer)!=kmers_in_ref.end())
		    {
		      //then the kmer does overlap a number, by definition.
		      first_del=false;//no need to ignore first SA interval (if it was in the nonvar bit would ignore)
		    }
		  else first_del=true;

		  bool precalc_done=true;
		  res_it = bidir_search_bwd(csa, (*it).first, (*it).second, 
					    (*it_rev).first, (*it_rev).second, 
					    p.begin(),p.begin()+p.size()-k, 
					    sa_intervals, sa_intervals_rev, 
					    sites, mask_a, maxx, first_del, precalc_done);

		  no_occ=0;

		  if (sa_intervals.size()<=10)
		    //proxy for mapping is "unique horizontally"
		    {
		      it=sa_intervals.begin();
		      no_occ=(*it).second-(*it).first; 
		      no_mapped++;
		      if (first_del==false) {
			assert(sites.front().empty());//becasue matches are all in non variable part of PRG
			repeats.clear();
			in_sites=0;
			for (auto ind=(*it).first;ind<(*it).second;ind++) {
			  if (mask_a[csa[ind]]!=0) {
			    in_sites++;
			    if (repeats.count(mask_s[csa[ind]])==0) repeats.insert(mask_s[csa[ind]]);
			    assert(mask_a[csa[ind]]==mask_a[csa[ind]+p.size()-1]);
			  }
			}
		      }
		     
		      it_ss=sites.begin();

		      while (it!=sa_intervals.end() && it_ss!=sites.end()) {
		       	if (it==sa_intervals.begin() && sites.front().empty()) {
			   assert(first_del==false);
			   for (auto ind=(*it).first;ind<(*it).second;ind++) 
			     if (mask_a[csa[ind]]!=0) {
			       covgs[(mask_s[csa[ind]]-5)/2][mask_a[csa[ind]]-1]=covgs[(mask_s[csa[ind]]-5)/2][mask_a[csa[ind]]-1]+1.0/(no_occ-in_sites+repeats.size()+sa_intervals.size()-1); //careful, might be dividing with more than we need to. size of sa_intervals is an overestimate of the number of horizontal matches, since a match that passed through 1st allele will be in a separate interval from other vertical matches from the same site
			       assert(mask_a[csa[ind]]==mask_a[csa[ind]+p.size()-1]);
			     }
			 }
			 else if ((it==sa_intervals.begin() && first_del==true) || (it!=sa_intervals.begin())) { //first_del=true - match in an interval starting with a number, all matches must be just to left of end marker
			  //if no_occ>1, two matches both starting at the end marker. If one crossed the start marker,
			  //sorina would have split into two SAs and here we are in one.
			  //so neither crosses the start marker, both start at the end. Since she only updates sites
			  //when you cross the left marker, it should be true that sites.front().back().second.size==0
			   auto it_s=*it_ss;
			   assert(!it_s.empty());
			   
			   invalid=false;
			   for (auto site_pair : it_s) {
			      auto site=site_pair.first;
			      auto allele=site_pair.second;
			      if (site_pair!=it_s.back() && allele.empty()) invalid=true; 
			   }
			   
			   if(!invalid) {
			     if (((*it).second-(*it).first)>1) assert(it_s.back().second.size()==0); //vertically non-unique 
			     for (auto site_pair : it_s) {
				auto site=site_pair.first;
				auto allele=site_pair.second;
				if (site_pair!=it_s.back() && site_pair!=it_s.front()) assert(allele.size()==1);
				if (site_pair==it_s.back() && allele.empty()) {
				  for (auto ind=(*it).first;ind<(*it).second;ind++) 
				    if (mask_a[csa[ind]]>0) {
				      if (first_del==false) {
					assert((no_occ-in_sites+repeats.size()+sa_intervals.size()-1)>0);
					covgs[(site-5)/2][mask_a[csa[ind]]-1]=covgs[(site-5)/2][mask_a[csa[ind]]-1]+1.0/(no_occ-in_sites+repeats.size()+sa_intervals.size()-1);
				      }
				      else
					covgs[(site-5)/2][mask_a[csa[ind]]-1]=covgs[(site-5)/2][mask_a[csa[ind]]-1]+1.0/sa_intervals.size();
				    }
				}
			        else if (!allele.empty()) {
				  assert(!allele.empty());
				  for (auto al:allele) {
				    if (first_del==false) {
				      assert((no_occ-in_sites+repeats.size()+sa_intervals.size()-1)>0);
				      covgs[(site-5)/2][al-1]=covgs[(site-5)/2][al-1]+1.0/(no_occ-in_sites+repeats.size()+sa_intervals.size()-1);
				    }
				    else
				      covgs[(site-5)/2][al-1]=covgs[(site-5)/2][al-1]+1.0/sa_intervals.size();
				  }
				}
			      }
			   }
			 }		       
		       ++it;
		       ++it_ss;
		      }
	          }
		  else 
		    {
		      no_occ=0;
		    }
		  sa_intervals.clear();
		  sa_intervals_rev.clear();
		  sites.clear();
	  }
	  else
	    {
	      no_occ=0;
	    }
	  //cout<<no_occ<<endl;
	  //clear p, sa_intervals etc
	  
	  no_reads++;
	  p.clear();
	}


	cout << "Finished mapping:"<<endl;
	timestamp();

	cout<<no_mapped<<endl;

	for (uint32_t i=0;i<covgs.size();i++) {
		for (uint32_t j=0;j<covgs[i].size();j++)
			out<<covgs[i][j]<<" ";
		out<<endl;
	}
	/*
	   for (int i=0;i<site_reads.size();i++) {
	   std::string name=string(_log)+"Site"+std::to_string(i);
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
