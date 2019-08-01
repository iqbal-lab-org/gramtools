#include <cctype>

#include "gtest/gtest.h"

#include "src_common/generate_prg.hpp"
#include "kmer_index/build.hpp"
#include "kmer_index/load.hpp"


using namespace gram;


TEST(DeserializeNextStats, GivenOneSearchStateWithThreePaths_CorrectlyIndexedKmerStats) {
    sdsl::int_vector<> kmers_stats = {3, 1, 42, 7};
    uint64_t stats_index = 0;
    auto result = deserialize_next_stats(stats_index, kmers_stats);
    IndexedKmerStats expected = {
            3,
            {1, 42, 7},
    };
    EXPECT_EQ(result, expected);
}


TEST(DeserializeNextStats, GivenTwoSearchStateWithMultiplePaths_CorrectlyIndexedKmerStats) {
    sdsl::int_vector<> kmers_stats = {3, 1, 42, 7, 2, 11, 33};
    uint64_t stats_index = 0;

    std::vector<IndexedKmerStats> all_stats = {};

    auto stats = deserialize_next_stats(stats_index, kmers_stats);
    all_stats.push_back(stats);
    stats_index += stats.count_search_states + 1;

    stats = deserialize_next_stats(stats_index, kmers_stats);
    all_stats.push_back(stats);

    const auto &result = all_stats;
    std::vector<IndexedKmerStats> expected = {
            {
                    3,
                    {1,  42, 7},
            },
            {
                    2,
                    {11, 33}
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(ParseSaIntervals, GivenOneKmerThreeSaIntervals_CorrectSearchStates) {
    KmerIndex kmer_index = {};
    sdsl::int_vector<3> all_kmers = {1, 2, 3, 4};
    sdsl::int_vector<> kmers_stats = {3, 0, 0, 0};

    Parameters parameters = {};
    parameters.sa_intervals_fpath = "@sa_intervals_fpath";
    parameters.kmers_size = 4;

    sdsl::int_vector<> sa_intervals = {42, 43, 52, 53, 62, 63};
    sdsl::util::bit_compress(sa_intervals);
    sdsl::store_to_file(sa_intervals, parameters.sa_intervals_fpath);

    parse_sa_intervals(kmer_index,
                       all_kmers,
                       kmers_stats,
                       parameters);
    const auto &result = kmer_index;
    KmerIndex expected = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {42, 43},
                            },
                            SearchState {
                                    SA_Interval {52, 53},
                            },
                            SearchState {
                                    SA_Interval {62, 63},
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(ParsePaths, GivenTwoPathsDifferentLengths_CorrectKmerIndex) {
    KmerIndex kmer_index = {};
    sdsl::int_vector<3> all_kmers = {1, 2, 3, 4};
    sdsl::int_vector<> kmers_stats = {2, 1, 2};

    Parameters parameters = {};
    parameters.paths_fpath = "@paths_fpath";
    parameters.kmers_size = 4;

    sdsl::int_vector<> paths = {42, 43, 52, 53, 62, 63};
    sdsl::util::bit_compress(paths);
    sdsl::store_to_file(paths, parameters.paths_fpath);

    parse_paths(kmer_index,
                all_kmers,
                kmers_stats,
                parameters);
    const auto &result = kmer_index;

    KmerIndex expected = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {},
                                    VariantSitePath {
                                            VariantLocus {42, 43}
                                    }
                            },
                            SearchState {
                                    SA_Interval {},
                                    VariantSitePath {
                                            VariantLocus {52, 53},
                                            VariantLocus {62, 63}
                                    }
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(ParseKmerIndex, GivenSingleKmerWithTwoSearchStates_CorrectKmerIndex) {
    Parameters parameters = {};
    parameters.kmers_size = 4;

    parameters.kmers_fpath = "@kmers_fpath";
    sdsl::int_vector<3> all_kmers = {1, 2, 3, 4};
    sdsl::store_to_file(all_kmers, parameters.kmers_fpath);

    parameters.kmers_stats_fpath = "@kmers_stats_fpath";
    sdsl::int_vector<> kmers_stats = {2, 1, 2};
    sdsl::util::bit_compress(kmers_stats);
    sdsl::store_to_file(kmers_stats, parameters.kmers_stats_fpath);

    parameters.sa_intervals_fpath = "@sa_intervals_fpath";
    sdsl::int_vector<> sa_intervals = {1, 1, 2, 2};
    sdsl::util::bit_compress(sa_intervals);
    sdsl::store_to_file(sa_intervals, parameters.sa_intervals_fpath);

    parameters.paths_fpath = "@paths_fpath";
    sdsl::int_vector<> paths = {42, 43, 52, 53, 62, 63};
    sdsl::util::bit_compress(paths);
    sdsl::store_to_file(paths, parameters.paths_fpath);

    auto result = kmer_index::load(parameters);

    KmerIndex expected = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {1, 1},
                                    VariantSitePath {
                                            VariantLocus {42, 43}
                                    }
                            },
                            SearchState {
                                    SA_Interval {2, 2},
                                    VariantSitePath {
                                            VariantLocus {52, 53},
                                            VariantLocus {62, 63}
                                    }
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}


TEST(ParseKmerIndex, GivenTwoKmersWithMultipleSearchStates_CorrectKmerIndex) {
    Parameters parameters = {};
    parameters.kmers_size = 4;

    parameters.kmers_fpath = "@kmers_fpath";
    sdsl::int_vector<3> all_kmers = {2, 2, 2, 2, 4, 4, 4, 4};
    sdsl::store_to_file(all_kmers, parameters.kmers_fpath);

    parameters.kmers_stats_fpath = "@kmers_stats_fpath";
    sdsl::int_vector<> kmers_stats = {1, 1, 2, 1, 2};
    sdsl::util::bit_compress(kmers_stats);
    sdsl::store_to_file(kmers_stats, parameters.kmers_stats_fpath);

    parameters.sa_intervals_fpath = "@sa_intervals_fpath";
    sdsl::int_vector<> sa_intervals = {1, 1, 1, 1, 2, 2};
    sdsl::util::bit_compress(sa_intervals);
    sdsl::store_to_file(sa_intervals, parameters.sa_intervals_fpath);

    parameters.paths_fpath = "@paths_fpath";
    sdsl::int_vector<> paths = {42, 43, 42, 43, 52, 53, 62, 63};
    sdsl::util::bit_compress(paths);
    sdsl::store_to_file(paths, parameters.paths_fpath);

    auto result = kmer_index::load(parameters);

    KmerIndex expected = {
            {{2, 2, 2, 2},
                    SearchStates {
                            SearchState {
                                    SA_Interval {1, 1},
                                    VariantSitePath {
                                            VariantLocus {42, 43}
                                    }
                            }
                    }
            },
            {{4, 4, 4, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {1, 1},
                                    VariantSitePath {
                                            VariantLocus {42, 43}
                                    }
                            },
                            SearchState {
                                    SA_Interval {2, 2},
                                    VariantSitePath {
                                            VariantLocus {52, 53},
                                            VariantLocus {62, 63}
                                    }
                            }
                    }
            }
    };
    EXPECT_EQ(result, expected);
}