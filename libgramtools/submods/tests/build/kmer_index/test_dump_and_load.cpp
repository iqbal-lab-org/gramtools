#include "gtest/gtest.h"

#include "src_common/common.hpp"
#include "build/kmer_index/build.hpp"
#include "build/kmer_index/dump.hpp"
#include "build/kmer_index/load.hpp"


using namespace gram;

/*********************/
/* Dumping (Writing) */
/*********************/

TEST(DumpKmers, GivenTwoKmers_CorrectAllKmersStructure) {
    BuildParams parameters = {};
    parameters.kmers_size = 4;
    parameters.kmers_fpath = "@kmers_fpath";

    KmerIndex kmer_index = {
            {{1, 2, 3, 4}, SearchStates {}},
            {{2, 4, 3, 4}, SearchStates {}},
    };

    sdsl::int_vector<3> all_kmers = dump_kmers(kmer_index, parameters);
    bool result = all_kmers == sdsl::int_vector<3> {1, 2, 3, 4, 2, 4, 3, 4}
                  or all_kmers == sdsl::int_vector<3> {2, 4, 3, 4, 1, 2, 3, 4};
    EXPECT_TRUE(result);
}

/***********/
/* Loading */
/***********/

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

/************************/
/* Dumping & loading */
/************************/

BuildParams setup_params(std::size_t const kmer_size){
    BuildParams parameters = {};
    parameters.kmers_size = kmer_size;
    parameters.kmers_fpath = "@kmers_fpath";
    parameters.kmers_stats_fpath = "@kmers_stats_fpath";
    parameters.sa_intervals_fpath = "@sa_intervals_fpath";
    parameters.paths_fpath = "@paths_fpath";

    return parameters;
}

TEST(DumpAndLoadIndex, SearchStatesWithNoVariants) {
    auto parameters = setup_params(4);

    KmerIndex kmer_index = {
            {{4, 4, 4, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {20000, 22000},
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
    ::kmer_index::dump(kmer_index, parameters);
    auto result = ::kmer_index::load(parameters);
    EXPECT_EQ(result, kmer_index);
}

TEST(DumpAndLoadIndex, SearchStateVariantsWithLargeIndices) {
    auto parameters = setup_params(4);

    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {1200000000, FIRST_ALLELE} // > 1 billion sites
                                    },
                                    VariantSitePath {},
                            },
                            SearchState {
                                    SA_Interval {7, 42},
                                    VariantSitePath {
                                            VariantLocus {5, 1200000000} // > 1 billion alleles
                                    },
                                    VariantSitePath {},
                            }
                    }
            }
    };

    ::kmer_index::dump(kmer_index, parameters);
    auto result = ::kmer_index::load(parameters);

    EXPECT_EQ(result, kmer_index);
}


TEST(DumpAndLoadIndex, TwoPathsWithMultipleElements) {
    auto parameters = setup_params(4);

    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, 1}
                                    },
                                    VariantSitePath {}
                            },
                            SearchState {
                                    SA_Interval {7, 42},
                                    VariantSitePath {
                                            VariantLocus {7, 3},
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {
                                        VariantLocus{9, ALLELE_UNKNOWN }
                                        }
                            }
                    }
            }
    };

    ::kmer_index::dump(kmer_index, parameters);
    auto result = ::kmer_index::load(parameters);

    EXPECT_EQ(result, kmer_index);
}


TEST(DumpAndLoadIndex, TwoKmersWithMultipleSearchStates) {
    auto parameters = setup_params(4);

    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, FIRST_ALLELE}
                                    },
                                    VariantSitePath {},
                            },
                            SearchState {
                                    SA_Interval {7, 7},
                                    VariantSitePath {
                                            VariantLocus {5, FIRST_ALLELE + 1}
                                    },
                                    VariantSitePath {},
                            },
                            SearchState {
                                    SA_Interval {8, 8},
                                    VariantSitePath {
                                            VariantLocus {5, FIRST_ALLELE + 1}
                                    },
                                    VariantSitePath {},
                            }
                    }
            },
            {{2, 4, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {9, 10},
                                    VariantSitePath {},
                                    VariantSitePath {},
                            },
                            SearchState {
                                    SA_Interval {11, 11},
                                    VariantSitePath {
                                            VariantLocus {5, FIRST_ALLELE + 1},
                                            VariantLocus {7, FIRST_ALLELE + 1}
                                    },
                                    VariantSitePath {},
                            }
                    }
            }
    };

    ::kmer_index::dump(kmer_index, parameters);
    auto result = ::kmer_index::load(parameters);

}

TEST(DumpAndLoadIndex, WithTraversingPaths) {
    auto parameters = setup_params(4);

    KmerIndex kmer_index = {
            { {1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, FIRST_ALLELE}
                                    },
                                    VariantSitePath{
                                            VariantLocus{7, ALLELE_UNKNOWN}
                                    },
                            },
                            SearchState {
                                    SA_Interval {7, 7},
                                    VariantSitePath {
                                            VariantLocus {5, FIRST_ALLELE + 1}
                                    },
                                    VariantSitePath {},
                            },
                            SearchState {
                                    SA_Interval {8, 8},
                                    VariantSitePath {
                                            VariantLocus {5, FIRST_ALLELE + 1}
                                    },
                                    VariantSitePath{
                                            VariantLocus{11, ALLELE_UNKNOWN},
                                            VariantLocus{9, ALLELE_UNKNOWN}
                                    },
                            }
                    }
            }
    };

    ::kmer_index::dump(kmer_index, parameters);
    auto result = ::kmer_index::load(parameters);
}
