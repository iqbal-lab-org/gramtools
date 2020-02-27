#include <cctype>

#include "gtest/gtest.h"

#include "src_common/common.hpp"
#include "kmer_index/build.hpp"
#include "kmer_index/dump.hpp"


using namespace gram;


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


TEST(DumpSaIntervals, GivenTwoSearchStates_CorrectSaIntervals) {
    BuildParams parameters = {};
    parameters.kmers_size = 4;
    parameters.kmers_fpath = "@kmers_fpath";
    parameters.sa_intervals_fpath = "@sa_intervals_fpath";

    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, 1}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {7, 42},
                                    VariantSitePath {
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };

    sdsl::int_vector<3> all_kmers = dump_kmers(kmer_index, parameters);
    auto stats = calculate_stats(kmer_index);
    dump_sa_intervals(stats, all_kmers, kmer_index, parameters);

    sdsl::int_vector<> result;
    sdsl::load_from_file(result, parameters.sa_intervals_fpath);

    sdsl::int_vector<> expected = {6, 6, 7, 42};
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(DumpPaths, GivenTwoPathsWithMultipleElements_CorrectSerializedPaths) {
    BuildParams parameters = {};
    parameters.kmers_size = 4;
    parameters.kmers_fpath = "@kmers_fpath";
    parameters.paths_fpath = "@paths_fpath";

    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, 1}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {7, 42},
                                    VariantSitePath {
                                            VariantLocus {7, 3},
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {
                                        VariantLocus{9, ALLELE_UNKNOWN }
                                        },
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };

    sdsl::int_vector<3> all_kmers = dump_kmers(kmer_index, parameters);
    auto stats = calculate_stats(kmer_index);
    dump_paths(stats, all_kmers, kmer_index, parameters);

    sdsl::int_vector<> result;
    sdsl::load_from_file(result, parameters.paths_fpath);

    sdsl::int_vector<> expected = {5, 1, 7, 3, 5, 2, 9, ALLELE_UNKNOWN};
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(DumpKmerEntryStats, GivenTwoKmersMultipleSearchStates_CorrectKmerEntryStats) {
    BuildParams parameters = {};
    parameters.kmers_size = 4;
    parameters.kmers_stats_fpath = "@kmers_stats_fpath";
    parameters.kmers_fpath = "@kmers_fpath";

    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, 1}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {7, 7},
                                    VariantSitePath {
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {8, 8},
                                    VariantSitePath {
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            },
            {{2, 4, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {9, 10},
                                    VariantSitePath {},
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {11, 11},
                                    VariantSitePath {
                                            VariantLocus {5, 2},
                                            VariantLocus {7, 2}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };

    sdsl::int_vector<3> all_kmers = dump_kmers(kmer_index, parameters);
    auto stats = calculate_stats(kmer_index);
    dump_kmers_stats(stats, all_kmers, kmer_index, parameters);

    sdsl::int_vector<> stats_kmer_entry;
    sdsl::load_from_file(stats_kmer_entry, parameters.kmers_stats_fpath);

    sdsl::int_vector<> expected_1 = {2, 0, 2, 3, 1, 1, 1};
    sdsl::util::bit_compress(expected_1);

    sdsl::int_vector<> expected_2 = {3, 1, 1, 1, 2, 0, 2};
    sdsl::util::bit_compress(expected_2);

    bool result = stats_kmer_entry == expected_1 or stats_kmer_entry == expected_2;
    EXPECT_TRUE(result);
}

TEST(DumpKmerEntryStats, GivenTraversingPaths_CorrectKmerEntryStats) {
    BuildParams parameters = {};
    parameters.kmers_size = 4;
    parameters.kmers_stats_fpath = "@kmers_stats_fpath";
    parameters.kmers_fpath = "@kmers_fpath";


    KmerIndex kmer_index = {
            { {1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantLocus {5, 1}
                                    },
                                    VariantSitePath{
                                            VariantLocus{7, ALLELE_UNKNOWN}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {7, 7},
                                    VariantSitePath {
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {8, 8},
                                    VariantSitePath {
                                            VariantLocus {5, 2}
                                    },
                                    VariantSitePath{
                                            VariantLocus{11, ALLELE_UNKNOWN},
                                            VariantLocus{9, ALLELE_UNKNOWN}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            }
    };

    sdsl::int_vector<3> all_kmers = dump_kmers(kmer_index, parameters);
    auto stats = calculate_stats(kmer_index);
    dump_kmers_stats(stats, all_kmers, kmer_index, parameters);

    sdsl::int_vector<> stats_kmer_entry;
    sdsl::load_from_file(stats_kmer_entry, parameters.kmers_stats_fpath);

    sdsl::int_vector<> expected = {3, 2, 1, 3};
    sdsl::util::bit_compress(expected);

    bool result = stats_kmer_entry == expected;
    EXPECT_TRUE(result);
}
