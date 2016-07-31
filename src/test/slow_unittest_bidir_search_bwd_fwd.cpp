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



TEST(BackwardSearchTest, NoVariantsSlowTest2){

  //PRG
  test_file2="../test_cases/36a.txt";

  //generate all substrings of PRG, use them all as queries
  ifstream ff(test_file2);
  ff >> prg;
  substrings=generate_all_substrings(prg);

  //dummy mask
  int a;
  mask_a.clear();
  for (a=0; a< prg.length(); a++)
    {
      mask_a.push_back(0); 
    }


  bool precalc=false;

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa 
    =csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  //  std::vector<int> mask_a;
 
  for (vector<string>::iterator it=substrings.begin();
       it<substrings.end();++it) {
    q=*it;
    
    bool first_del=false;
    bool precalc = false;
    int occ_expt=0;
    int pos=prg.find(q,0);

    while (pos!=string::npos) {
      occ_expt++;
      pos=prg.find(q,pos+1);
    }

    for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
    }

    std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del, precalc);

    uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_TRUE(first_del==false);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();

    csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
    first_del=false;
    res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del, precalc);  

    no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_TRUE(first_del==false);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();
    p.clear();
  }
}


TEST(BackwardSearchTest, NoVariantsABCABCTest3){

  //PRG
  test_file2="../test_cases/abc_abc_abc.txt";

  //generate all substrings of PRG, use them all as queries
  ifstream ff(test_file2);
  ff >> prg;
  substrings=generate_all_substrings(prg);

  //dummy mask
  int a;
  mask_a.clear();
  for (a=0; a< prg.length(); a++)
    {
      mask_a.push_back(0); 
    }


  bool precalc=false;

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa 
    =csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  //  std::vector<int> mask_a;
 
  for (vector<string>::iterator it=substrings.begin();it<substrings.end();++it) {
    q=*it;

    bool first_del=false;
    bool precalc = false;
    int occ_expt=0;
    int pos=prg.find(q,0);

    while (pos!=string::npos) {
      occ_expt++;
      pos=prg.find(q,pos+1);
    }

    for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
    }

    std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del, precalc);

    uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_TRUE(first_del==false);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();

    csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
    first_del=false;
    res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del, precalc);  

    no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_TRUE(first_del==false);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();
    p.clear();
  }
}


TEST(BackwardSearchTest, NoVariantsACTG4){

  //PRG
  test_file2="../test_cases/actg.txt";

  //generate all substrings of PRG, use them all as queries
  ifstream ff(test_file2);
  ff >> prg;
  substrings=generate_all_substrings(prg);

  //dummy mask
  int a;
  mask_a.clear();
  for (a=0; a< prg.length(); a++)
    {
      mask_a.push_back(0); 
    }


  bool precalc=false;

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa 
    =csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  //  std::vector<int> mask_a;
 
  for (vector<string>::iterator it=substrings.begin();it<substrings.end();++it) {
    q=*it;

    bool first_del=false;
    bool precalc = false;
    int occ_expt=0;
    int pos=prg.find(q,0);

    while (pos!=string::npos) {
      occ_expt++;
      pos=prg.find(q,pos+1);
    }

    for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
    }

    std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del, precalc);

    uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_TRUE(first_del==false);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();

    csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
    first_del=false;
    res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del, precalc);  

    no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_TRUE(first_del==false);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();
    p.clear();
  }
}


TEST(BackwardSearchTest, NoVariantsSlow_MSP34_200bp_Test5){

  //PRG
  test_file2="../test_cases/MSP3.4_200_bases.txt";

  //generate all substrings of PRG, use them all as queries
  ifstream ff(test_file2);
  ff >> prg;
  substrings=generate_all_substrings(prg);

  //dummy mask
  int a;
  mask_a.clear();
  for (a=0; a< prg.length(); a++)
    {
      mask_a.push_back(0); 
    }


  bool precalc=false;

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa 
    =csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",true, false);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  //  std::vector<int> mask_a;
 
  for (vector<string>::iterator it=substrings.begin();it<substrings.end();++it) {
    q=*it;

    bool first_del=false;
    bool precalc = false;
    int occ_expt=0;
    int pos=prg.find(q,0);

    while (pos!=string::npos) {
      occ_expt++;
      pos=prg.find(q,pos+1);
    }

    for (uint16_t i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
    }

    std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size(),0,csa.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del, precalc);

    uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_TRUE(first_del==false);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();

    csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> csa_rev=csa_constr(test_file2, "int_alphabet_file","memory_log_file","csa_file",false, false);
    first_del=false;
    res_it=bidir_search_fwd(csa_rev,0,csa_rev.size(),0,csa_rev.size(),p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del, precalc);  

    no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
    EXPECT_TRUE(first_del==false);
    EXPECT_EQ(1,sa_intervals.size());
    EXPECT_EQ(no_occ,occ_expt);

    sa_intervals.clear();
    sa_intervals_rev.clear();
    sites.clear();
    p.clear();
  }
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
