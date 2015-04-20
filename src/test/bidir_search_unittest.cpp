#include "gtest/gtest.h"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <iostream>
#include <vector>
#include <string>

using namespace sdsl;
using namespace std;

class BidirSearchTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    construct(csa,test_file,4);
    construct(csa_rev,test_file_rev,4);
    l=0;
    r=csa.size()-1;
  }
  virtual void TearDown() {
    util::clear(csa);
  }

  csa_wt<wt_int<rrr_vector<63>>> csa, csa_rev;
  uint32_t l, r,l_rev,r_rev;

};

TEST_F(BidirSearchTest, fwd_ind){
  l_old=l
  r_old=r
  occ=0 
  ASSERT_TRUE("c==1 or c==2 or c==3 or c==4")

  bidir_search(csa,l,r,l_rev,r_rev,c);

  EXPECT_EQ(r_rev-l_rev, r-l);
  ASSERT_LT(l,r);
  ASSERT_GE(l,csa.C[csa.char2comp[c]]);
  if (c<max) ASSERT_LE(r,csa.C[csa.char2comp[(c+1)]]);

  for (int i=l_old;i<r_old;i++)
    if (text[csa[i]-1]==c) {
      occ++;
      rnk=l+occ-1-csa.C[csa.char2comp[c]];
      EXPECT_EQ(i,csa.bwt.select(rnk));
    }

  EXPECT_EQ(occ,r-l);
  EXPECT_EQ(occ,r_rev-l_rev);
}

TEST_F(BidirSearchTest, bwd_ind){
  l_rev_old=l_rev
  r_rev_old=r_rev
  ASSERT_TRUE("c==1 or c==2 or c==3 or c==4")

  bidir_search(csa,l,r,l_rev,r_rev,c)
  
  ASSERT_LT(l_rev,r_rev);
  ASSERT_LE(l_rev_old,l_rev);
  ASSERT_GE(r_rev_old,r_rev);


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  if (argc!=2) {
    cout << "Usage: " << argv[0] << " test_file" << endl;
    cout << " (1) Reverses test_file; stores it in test_file_rev." << endl;
    cout << " (2) Performs tests." << endl;
    cout << " (3) Deletes test_file_reverse." << endl;
    return 1;
  }

  test_file = argv[1];
  test_file_rev = test_file + "_rev";

  int_vector<4> text;
  load_vector_from_file(text, test_file, 1);
  size_type n = text.size();
  int_vector<4> text_rev(n);
  uint32_t max=0;
  for (size_type i=0; i<n; i++) {
    text_rev[n-1-i] = text[i];
    if (text[i]>max) max=text[i];
  }
  char* text2 = (char*)text_rev.data();
  ofstream of(test_file_rev, ofstream::binary);
  of.write(text2, n);
  of.close();

  return RUN_ALL_TESTS();
}

