#include <sstream>

#include <boost/iostreams/copy.hpp>         // copy
#include <boost/iostreams/filter/gzip.hpp>  // decompressor
#include <boost/iostreams/filtering_streambuf.hpp>
#include "gtest/gtest.h"

#include "build/check_ref.hpp"
#include "prg/linearised_prg.hpp"

using namespace gram;

coverage_Graph setup_cov_graph(std::string const& prg_string) {
  auto markers = prg_string_to_ints(prg_string);
  PRG_String p{markers};
  return coverage_Graph{p};
}

TEST(GetFirstPrgPath, NonNestedPrg) {
  auto cov_graph = setup_cov_graph("[AC,GG]GG[A,T,C]CA[,G]C");
  auto ref_path = PrgRefChecker::get_first_prg_path(cov_graph);
  EXPECT_EQ(ref_path, std::string{"ACGGACAC"});
}

TEST(GetFirstPrgPath, NonNestedPrg2) {
  auto cov_graph = setup_cov_graph("A[AAA,GG]GG[A,]CAC");
  auto ref_path = PrgRefChecker::get_first_prg_path(cov_graph);
  EXPECT_EQ(ref_path, std::string{"AAAAGGACAC"});
}

TEST(GetFirstPrgPath, NestedPrg) {
  auto cov_graph = setup_cov_graph("[AC[CG,C]TTT[C[A,G],G]T,GG]CA[A,G[A,C]]C");
  auto ref_path = PrgRefChecker::get_first_prg_path(cov_graph);
  EXPECT_EQ(ref_path, std::string{"ACCGTTTCATCAAC"});
}

class TestRefMatchesFirstPrgPath : public ::testing::Test {
 protected:
  void SetUp() {
    // First path: "AACTCCAAACG"
    cov_graph = setup_cov_graph("A[AC,TT]TCC[AAA[C,A],G]G");
  }
  coverage_Graph cov_graph;
};

TEST_F(TestRefMatchesFirstPrgPath, CorrectRef_Passes) {
  std::istringstream ss{"AACTCCAAACG"};
  PrgRefChecker(ss, cov_graph);
}

TEST_F(TestRefMatchesFirstPrgPath,
       CorrectRefWithFastaHeader_HeaderIgnoredAndPasses) {
  std::istringstream ss{">chrom1\nAACTCCAAACG"};
  PrgRefChecker(ss, cov_graph);
}

TEST_F(TestRefMatchesFirstPrgPath, IncorrectRef_Fails) {
  std::istringstream ss{"ATTTTTTT"};
  EXPECT_THROW(PrgRefChecker(ss, cov_graph), std::runtime_error);
}

TEST_F(TestRefMatchesFirstPrgPath, LowerCaseCorrectRef_Passes) {
  std::istringstream ss{"aactccaaacg"};
  PrgRefChecker(ss, cov_graph);
}

TEST_F(TestRefMatchesFirstPrgPath, GzipCompressedCorrectRef_Passes) {
  using namespace boost::iostreams;
  std::stringstream ss{"AACTCCAAACG"};
  filtering_istreambuf compressor;
  compressor.push(gzip_compressor());
  compressor.push(ss);

  std::stringstream compressed;
  boost::iostreams::copy(compressor, compressed);
  PrgRefChecker(compressed, cov_graph, true);  // Say input stream is gzipped
}
