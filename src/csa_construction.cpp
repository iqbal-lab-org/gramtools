#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <vector>
#include <cstdint>
#include <iostream>
#include <fstream>


using namespace sdsl;
using namespace std;

 csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_constr(std::string fname, std::vector<std::vector<int>>& covgs, char* int_al_fname) {
   std::ifstream f(fname);
   std::string prg;
   std::ofstream out("csa-mem.html");
   std::streambuf *coutbuf = std::cout.rdbuf();
   FILE* fp;

   uint64_t *prg_int=(uint64_t*)malloc(prg.length()*sizeof(uint64_t));
  //always check if the malloc has succeeded
   if (prg_int==NULL)
     {
      //die
     }
   int i=0;
   int ii=0;
   while (i<prg.length()) {
       if (isdigit(prg[i])) {
	 int j=1;
	 while(isdigit(prg[i+1])) {
	   j++;
	   i++;
	 }
	 auto al_ind=prg.substr(i-j+1,j);
	 //uint64_t l=(uint64_t) stoull(al_ind,NULL,0);
	 auto l=stoull(al_ind,NULL,0);
	 //uint64_t l=boost::lexical_cast<uint64_t>(al_ind); 
	 prg_int[ii]=l;
       }
       else {
	 if (prg[i]=='A') prg_int[ii]=1;
	 if (prg[i]=='C') prg_int[ii]=2;
	 if (prg[i]=='G') prg_int[ii]=3;
	 if (prg[i]=='T') prg_int[ii]=4;
       }
       i++;
       ii++;
   }

   fp=fopen(int_al_fname,"wb");
   fwrite(prg_int,sizeof(uint64_t),ii,fp);
   fclose(fp);
   std::cout.rdbuf(out.rdbuf());   
 
   memory_monitor::start();
   csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa;
   construct(csa,int_al_fname,8);  
   memory_monitor::stop();
   memory_monitor::write_memory_log<HTML_FORMAT>(cout);
 
   std::cout.rdbuf(coutbuf);
  //  store_to_file(csa,"CSA_file");
   return(csa);
}
