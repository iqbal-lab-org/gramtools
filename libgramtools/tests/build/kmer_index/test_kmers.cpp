#include "build/kmer_index/kmers.hpp"
#include "gtest/gtest.h"
#include "submod_resources.hpp"

using namespace gram::submods;

TEST(GetReversedKmers,
     GivenRandomlyArrangedReverseKmers_KmersReversedAndSortedByRightMostBase) {
  ordered_vector_set<Sequence> kmers = {
      {2, 4, 1},
      {1, 3, 5},
      {1, 3, 4},
      {3, 4, 5},
  };

  std::vector<Sequence> result = reverse(kmers);
  Sequences expected = {
      {4, 3, 1},
      {5, 3, 1},
      {1, 4, 2},
      {5, 4, 3},
  };
  EXPECT_EQ(result, expected);
}

TEST(GetReversedKmers, GivenSingleReverseKmer_CorrectReversedKmer) {
  ordered_vector_set<Sequence> kmers = {
      {2, 4, 1},
  };

  std::vector<Sequence> result = reverse(kmers);
  Sequences expected = {
      {1, 4, 2},
  };
  EXPECT_EQ(result, expected);
}

TEST(GetReversedKmers, SortingReverseKmerFromRightToLeft_CorrectReversedKmers) {
  ordered_vector_set<Sequence> kmers = {
      {1, 3, 5},
      {2, 4, 1},
  };

  std::vector<Sequence> result = reverse(kmers);
  Sequences expected = {
      {5, 3, 1},
      {1, 4, 2},
  };
  EXPECT_EQ(result, expected);
}

TEST(GetPrefixDiffs, GivenKmersDifferInLeftMostBaseOnly_CorrectPrefixDiffs) {
  std::vector<Sequence> kmers = {
      {1, 3, 1},
      {2, 3, 1},
      {3, 3, 1},
      {4, 3, 1},
  };

  auto result = get_prefix_diffs(kmers);
  std::vector<Sequence> expected = {
      {1, 3, 1},
      {2},
      {3},
      {4},
  };
  EXPECT_EQ(result, expected);
}

TEST(GetPrefixDiffs, GivenKmerDifferInRightMostBaseOnly_CorrectPrefixDiffs) {
  std::vector<Sequence> kmers = {
      {1, 3, 1},
      {2, 3, 1},
      {1, 3, 2},
  };

  auto result = get_prefix_diffs(kmers);
  std::vector<Sequence> expected = {
      {1, 3, 1},
      {2},
      {1, 3, 2},
  };
  EXPECT_EQ(result, expected);
}

TEST(GetPrefixDiffs, GivenMixOfOrderedKmers_CorrectPrefixDiffs) {
  std::vector<Sequence> kmers = {
      {1, 3, 1}, {2, 3, 1}, {1, 3, 2}, {1, 4, 2}, {3, 4, 2},
  };

  auto result = get_prefix_diffs(kmers);
  std::vector<Sequence> expected = {
      {1, 3, 1}, {2}, {1, 3, 2}, {1, 4}, {3},
  };
  EXPECT_EQ(result, expected);
}

TEST(GetAllKmers, GenerateAllKmersLengthThree_CorrectOrder) {
  auto result = get_all_kmers(3);

  std::vector<Sequence> expected = {
      {1, 1, 1}, {2, 1, 1}, {3, 1, 1}, {4, 1, 1}, {1, 2, 1}, {2, 2, 1},
      {3, 2, 1}, {4, 2, 1}, {1, 3, 1}, {2, 3, 1}, {3, 3, 1}, {4, 3, 1},
      {1, 4, 1}, {2, 4, 1}, {3, 4, 1}, {4, 4, 1}, {1, 1, 2}, {2, 1, 2},
      {3, 1, 2}, {4, 1, 2}, {1, 2, 2}, {2, 2, 2}, {3, 2, 2}, {4, 2, 2},
      {1, 3, 2}, {2, 3, 2}, {3, 3, 2}, {4, 3, 2}, {1, 4, 2}, {2, 4, 2},
      {3, 4, 2}, {4, 4, 2}, {1, 1, 3}, {2, 1, 3}, {3, 1, 3}, {4, 1, 3},
      {1, 2, 3}, {2, 2, 3}, {3, 2, 3}, {4, 2, 3}, {1, 3, 3}, {2, 3, 3},
      {3, 3, 3}, {4, 3, 3}, {1, 4, 3}, {2, 4, 3}, {3, 4, 3}, {4, 4, 3},
      {1, 1, 4}, {2, 1, 4}, {3, 1, 4}, {4, 1, 4}, {1, 2, 4}, {2, 2, 4},
      {3, 2, 4}, {4, 2, 4}, {1, 3, 4}, {2, 3, 4}, {3, 3, 4}, {4, 3, 4},
      {1, 4, 4}, {2, 4, 4}, {3, 4, 4}, {4, 4, 4},
  };

  for (uint64_t i = 0; i < result.size(); ++i) {
    auto result_kmer = result[i];
    auto expected_kmer = expected[i];
    EXPECT_EQ(result_kmer, expected_kmer);
  }
}

TEST(GenerateKmers, GenerateAllKmersOfSizeThree_CorrectSpotCheck) {
  auto kmers = generate_all_kmers(3);
  std::vector<Sequence> expected_kmers = {
      {1, 1, 1}, {1, 1, 2}, {1, 1, 3}, {1, 1, 4}, {1, 2, 1},
      {1, 2, 2}, {1, 2, 3}, {1, 2, 4}, {1, 3, 1}, {3, 3, 3},
      {4, 4, 2}, {1, 4, 2}, {4, 4, 4},
  };

  for (const auto &expected_kmer : expected_kmers) {
    auto result =
        std::find(kmers.begin(), kmers.end(), expected_kmer) != kmers.end();
    EXPECT_TRUE(result);
  }
}
