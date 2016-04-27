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

string test_file,q,test_file2,query,mask_file;
std::vector<uint8_t> p;
std::vector<std::vector<int>> covgs;
string prg,prg2;
vector<string> substrings;
std::vector<int> mask_a;

/*
TEST(BackwardSearchTest, NoVariants){

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  std::vector<int> mask_a;
 
  for (vector<string>::iterator it=substrings.begin();it<substrings.end();++it) {
    q=*it;

    bool first_del=false;
    int occ_expt=0;
    int pos=prg.find(q,0);

    while (pos!=string::npos) {
      occ_expt++;
      pos=prg.find(q,pos+1);
    }

    for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
    }

    std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del);

    uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_EQ(false,first_del);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();

    csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
    first_del=false;
    res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del);  

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

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,6,first_del);

  uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,6,first_del);  

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

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);

  uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);  

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

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,6,first_del);

  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
    no_occ+=(*it).second-(*it).first;
  
  EXPECT_EQ(false,first_del);
  EXPECT_EQ(2,sa_intervals.size());
  EXPECT_EQ(no_occ,2);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,6,first_del);  

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

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);

  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);  

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

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);

  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(false,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,2);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);  

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

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);
  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);  

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

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);
  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(false,first_del);
  EXPECT_EQ(3,sa_intervals.size());
  EXPECT_EQ(no_occ,3);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file2,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,8,first_del);  

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
*/

TEST(BackwardSearchTest, One_match_many_sites){

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::pair<uint64_t,uint64_t>>::iterator it;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  bool first_del=false;

  q=query;
  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,16,first_del);
  uint64_t no_occ=0;
  for (it=sa_intervals.begin();it!=sa_intervals.end();++it) 
    no_occ+=(*it).second-(*it).first;

  EXPECT_EQ(true,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,1);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file2,"int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,16,first_del);  

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

  test_file = argv[1];
  // q=argv[2];
  ifstream f(test_file);
  f>>prg;
  substrings=generate_all_substrings(prg);

  test_file2=argv[2];
  query=argv[3];
  
  mask_file=argv[4];
  ifstream g(mask_file);

  int a;
  while (g>>a) mask_a.push_back(a);

  return RUN_ALL_TESTS();
}
