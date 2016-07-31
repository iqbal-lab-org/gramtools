#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "gtest/gtest.h"
#include "bwt_search.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace sdsl;
using namespace std;

string q,test_file2,query,mask_file;
std::vector<uint8_t> p;
std::vector<std::vector<int>> covgs;
string prg,prg2;
vector<string> substrings;
std::vector<int> mask_a;

//forward declare
vector<string> generate_all_substrings(string q);

TEST(BackwardSearchTest, NoVariants1){

  //PRG
  test_file2="../test_cases/one_byte.txt";

  //generate all substrings of PRG, use them all as queries
  string temp;
  ifstream ff(test_file2);
  ff >> temp;//in this case, is a one char PRG, no need for substrings
  substrings.push_back(temp);

  //dummy mask
  int a;
  mask_a.clear();
  for (a=0; a< temp.length(); a++)
    {
      mask_a.push_back(0);
    }


  bool precalc=false;

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  //  std::vector<int> mask_a;
 
  for (vector<string>::iterator it=substrings.begin();it<substrings.end();++it) {
    q=*it;

    bool first_del=false;
    bool precalc = false;
    int occ_expt=1;

    for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
    }

    std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,4,first_del, precalc);

    uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_EQ(false,first_del);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();

    csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
    first_del=false;
    res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,4,first_del, precalc);  

    no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_EQ(false,first_del);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();
    p.clear();
  }

}







TEST(BackwardSearchTest, OneSNP){


  test_file2="../test_cases/one_snp.txt";
  query="ttacacagaactagagag";
  mask_file="../test_cases/one_snp_mask_a.txt";
  ifstream g(mask_file);
  bool precalc=false;

  int a;
  mask_a.clear();
  while (g>>a) mask_a.push_back(a);

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=csa_constr(test_file2,"int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,6,first_del, precalc);

  uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,6,first_del,precalc);  

  no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();
  p.clear();
}


TEST(BackwardSearchTest, TwoSNPs){

  test_file2="../test_cases/two_snps.txt";
  query="ttacacagaactagaagcag";
  mask_file="../test_cases/two_snps_mask_a.txt";
  ifstream g(mask_file);

  int a;
  mask_a.clear();
  while (g>>a) mask_a.push_back(a);

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;
  bool precalc=false;
  q=query;
  for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),
							 p.begin(),p.end(), 
							 sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);

  uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);  

  no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();
  p.clear();
}


TEST(BackwardSearchTest, Two_matches_one_variable_one_nonvariable_region){

  test_file2="../test_cases/two_matches_var_nonvar.txt";
  query="acagaac";
  mask_file="../test_cases/two_matches_var_nonvar_mask_a.txt";
  ifstream g(mask_file);

  int a;
  mask_a.clear();
  while (g>>a) mask_a.push_back(a);

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;
  bool precalc=false;
  q=query;
  for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,6,first_del, precalc);

  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
    no_occ+=(*it).second-(*it).first;
  
  EXPECT_EQ(false,first_del);
  EXPECT_EQ(2,sa_intervals.size());
  EXPECT_EQ(no_occ,2);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,6,first_del, precalc);  

  no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(false,first_del);
  EXPECT_EQ(2,sa_intervals.size());
  EXPECT_EQ(no_occ,2);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();
  p.clear();
}


TEST(BackwardSearchTest, Two_long_sites){


  test_file2="../test_cases/two_long_sites.txt";
  query="gctcggctcgatgactagatagatagcgaggcaac";
  mask_file="../test_cases/two_long_sites_mask_a.txt";
  ifstream g(mask_file);

  int a;
  mask_a.clear();
  while (g>>a) mask_a.push_back(a);

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=
    csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;
  bool precalc=false;
  q=query;
  for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);

  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);  

  no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();
  p.clear();
}


TEST(BackwardSearchTest, Match_within_long_site_match_outside){



  test_file2="../test_cases/match_within_long_site.txt";
  //  query="tagacacacagtgtcgcctcgtcggctttgagtggtgctagacctca";
  query="ctgctccacacagaga";
  mask_file="../test_cases/match_within_long_site_mask_a.txt";
  ifstream g(mask_file);

  int a;
  mask_a.clear();
  while (g>>a) mask_a.push_back(a);

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=
    csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;
  bool precalc=false;
  q=query;
  for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);

  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(false,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,2);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
  first_del=false;

  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);  

  no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(false,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,2);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();
  p.clear();
}



TEST(BackwardSearchTest, Long_site_and_repeated_snp_on_edge_of_site){

  test_file2="../test_cases/repeated_snp_on_both_edges.txt";
  query="tagacacacagtgtcgcctcgtcggctttgagtggtgctagacccca";
  mask_file="../test_cases/match_within_long_site_mask_a.txt";
  ifstream g(mask_file);

  int a;
  mask_a.clear();
  while (g>>a) mask_a.push_back(a);


  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true,false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;
  bool precalc=false;

  q=query;
  for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);
  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false,false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);  

  no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();
  p.clear();
}




TEST(BackwardSearchTest, Multiple_matches_over_multiple_sites){

  test_file2="../test_cases/multiple_matches_multiple_sites.txt";
  query="tgata";
  mask_file="../test_cases/multiple_matches_multiple_sites_mask_a.txt";
  ifstream g(mask_file);

  int a;
  mask_a.clear();
  while (g>>a) mask_a.push_back(a);


  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;
  bool precalc=false;
  q=query;
  for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);
  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(false,first_del);
  EXPECT_EQ(3,sa_intervals.size());
  EXPECT_EQ(no_occ,3);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del, precalc);  

  no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(false,first_del);
  EXPECT_EQ(3,sa_intervals.size());
  EXPECT_EQ(no_occ,3);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();
  p.clear();
}

TEST(BackwardSearchTest, One_match_many_sites){


  test_file2="../test_cases/One_match_many_sites.txt";
  query="cctacacatgatcgtgatcaccatagaggtcgctgggtccat";
  mask_file="../test_cases/One_match_many_sites_mask_a.txt";
  ifstream g(mask_file);

  int a;
  mask_a.clear();
  while (g>>a) mask_a.push_back(a);

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;
  bool precalc=false;
  q=query;
  for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,16,first_del, precalc);
  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2,"int_alphabet_file","memory_log_file","csa_file",false, false);
  first_del=false;

  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,16,first_del, precalc);  

  no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();
  p.clear();
}



vector<string> generate_all_substrings(string q) {

  vector<string> substrings;
  int n,c,i;

  n=q.size();
  for (c=0; c<n; c++) {
    for (i=1; i<=n-c; i++) {
      substrings.push_back(q.substr(c,c+i));
    }
  }
  return(substrings);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);


  return RUN_ALL_TESTS();
}
