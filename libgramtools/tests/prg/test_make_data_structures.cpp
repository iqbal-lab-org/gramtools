#include "gtest/gtest.h"

#include "prg/make_data_structures.hpp"
#include "submod_resources.hpp"

using namespace gram::submods;

TEST(GetNumVarSites, NoSites) {
  auto prg_raw = encode_prg("c");
  auto prg_info = generate_prg_info(prg_raw);
  auto result = prg_info.num_variant_sites;
  uint64_t expected = 0;
  EXPECT_EQ(result, expected);
}

TEST(GetNumVarSites, UnNestedPrgString) {
  auto prg_raw = encode_prg("a5g6t6cccc11g12tttt12");
  auto prg_info = generate_prg_info(prg_raw);
  auto result = prg_info.num_variant_sites;
  auto expected = 2;
  EXPECT_EQ(result, expected);
}

TEST(GetNumVarSites, Nested_PrgString) {
  auto prg_raw = prg_string_to_ints("[[A,C,G]A,T]T[,C][GA,CT]");
  auto prg_info = generate_prg_info(prg_raw);
  auto result = prg_info.num_variant_sites;
  uint64_t expected = 4;
  EXPECT_EQ(result, expected);
}

TEST(BuildChildMap, GivenParentalMap_CorrectChildMap) {
  // Site 5 has two sites nested in haplogroup 1, and one in haplogroup 2.
  // Note: parental_map / quasimap stores allele haplogroups as 1-based,
  // but child_map moves them to 0-based (consistent with infer).
  parental_map par_map{
      {7, VariantLocus{5, FIRST_ALLELE}},
      {9, VariantLocus{5, FIRST_ALLELE}},
      {11, VariantLocus{5, FIRST_ALLELE + 1}},
      {15, VariantLocus{13, FIRST_ALLELE + 2}},
  };

  auto result = build_child_map(par_map);
  // Sort the internal vectors to be independent from parental hash map element
  // ordering
  for (auto& entry : result) {
    for (auto& entry2 : entry.second) {
      std::sort(entry2.second.begin(), entry2.second.end());
    }
  }
  child_map expected{{5,
                      haplo_map{
                          {0, marker_vec{7, 9}},
                          {1, marker_vec{11}},
                      }},
                     {13, haplo_map{{2, marker_vec{15}}}}};

  EXPECT_EQ(result, expected);
}
