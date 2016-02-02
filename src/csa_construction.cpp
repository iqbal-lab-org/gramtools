#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <vector>
#include <cstdint>
#include <iostream>
#include <fstream>


using namespace sdsl;
using namespace std;

//make SA sampling density and ISA sampling density customizable
// what is this fname - the binary file we create of the CSA?
csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_constr(std::string fname, 
   std::vector<std::vector<int>>& covgs, 
   char* int_al_fname, 
   char* memory_log_fname, 
   char* csa_file) 

{
   std::ifstream f(fname);
   std::string prg;
   std::ofstream out(memory_log_fname);
   std::streambuf *coutbuf = std::cout.rdbuf();
   FILE* fp;

   f>>prg;
   f.close();

   uint64_t *prg_int=(uint64_t*)malloc(prg.length()*sizeof(uint64_t));

 //always check if the malloc has succeeded
   if (prg_int==NULL)
     {
       exit(1);
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
       else 
       {
        	 if (prg[i]=='A' or prg[i]=='a') prg_int[ii]=1;
	         if (prg[i]=='C' or prg[i]=='c') prg_int[ii]=2;
	         if (prg[i]=='G' or prg[i]=='g') prg_int[ii]=3;
	         if (prg[i]=='T' or prg[i]=='t') prg_int[ii]=4;
        }
       i++;
       ii++;// so ii keeps track of actual base position - it's aware of numbers with more than one digit
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
   store_to_file(csa,csa_file);
   return(csa);
}
