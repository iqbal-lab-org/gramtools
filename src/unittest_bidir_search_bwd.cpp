#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "gtest/gtest.h"
#include "bwt_search.h"
#include <iostream>
#include <vector>
#include <string>

using namespace sdsl;
using namespace std;

string test_file,q;
ifstream f(test_file);
std::vector<uint8_t> p;
std::vector<std::vector<int>> covgs;


TEST(BackwardSearchTest, NoVariants){

  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(test_file,covgs, "int_alphabet_file","memory_log_file","csa_file",true);

  std::list<std::pair<uint64_t,uint64_t>> sa_intervals, sa_intervals_rev;
  std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
  std::vector<int> mask_a;
  bool first_del=false;

  for (int i=0;i<q.length();i++) {
       if (q[i]=='A' or q[i]=='a') p.push_back(1);
       if (q[i]=='C' or q[i]=='c') p.push_back(2);
       if (q[i]=='G' or q[i]=='g') p.push_back(3);
       if (q[i]=='T' or q[i]=='t') p.push_back(4);
  }

  
  std::vector<uint8_t>::iterator res_it=bidir_search_bwd(csa,0,csa.size()-1,0,csa.size()-1,p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del);

  uint64_t no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
  cout<<(*sa_intervals.begin()).second<<" "<<csa[(*sa_intervals.begin()).second-1]<<endl;
  cout<<csa[2]<<endl;
  cout<<(*sa_intervals.begin()).first<<" "<<csa[(*sa_intervals.begin()).first]<<endl;

  EXPECT_EQ(false,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,3);

  sa_intervals.clear();
  sa_intervals_rev.clear();
  sites.clear();


  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_rev=csa_constr(test_file,covgs, "int_alphabet_file","memory_log_file","csa_file",false);
  first_del=false;
  res_it=bidir_search_fwd(csa_rev,0,csa_rev.size()-1,0,csa_rev.size()-1,p.begin(),p.end(), sa_intervals,sa_intervals_rev,sites,mask_a,5,first_del);  

  cout<<(*sa_intervals.begin()).second<<" "<<csa_rev[(*sa_intervals.begin()).second-1]<<endl;
  cout<<(*sa_intervals.begin()).first<<" "<<csa_rev[(*sa_intervals.begin()).first]<<endl;

  no_occ=(*sa_intervals.begin()).second-(*sa_intervals.begin()).first;
  EXPECT_EQ(false,first_del);
  EXPECT_EQ(1,sa_intervals.size());
  EXPECT_EQ(no_occ,3);

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  test_file = argv[1];
  q=argv[2];

  return RUN_ALL_TESTS();
}
