#include <cctype>

#include "gtest/gtest.h"

#include "../test_utils.hpp"
#include "kmer_index/build.hpp"
#include "kmer_index/dump.hpp"


using namespace gram;


TEST(DumpKmers, GivenTwoKmers_CorrectAllKmersStructure) {
    Parameters parameters = {};
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
    Parameters parameters = {};
    parameters.kmers_size = 4;
    parameters.kmers_fpath = "@kmers_fpath";
    parameters.sa_intervals_fpath = "@sa_intervals_fpath";

    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantSite {5, 1}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {7, 42},
                                    VariantSitePath {
                                            VariantSite {5, 2}
                                    },
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
    Parameters parameters = {};
    parameters.kmers_size = 4;
    parameters.kmers_fpath = "@kmers_fpath";
    parameters.paths_fpath = "@paths_fpath";

    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantSite {5, 1}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {7, 42},
                                    VariantSitePath {
                                            VariantSite {5, 2},
                                            VariantSite {7, 3}
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

    sdsl::int_vector<> expected = {5, 1, 5, 2, 7, 3};
    sdsl::util::bit_compress(expected);
    EXPECT_EQ(result, expected);
}


TEST(DumpKmerEntryStats, GivenTwoKmersMultipleSearchStates_CorrectKmerEntryStats) {
    Parameters parameters = {};
    parameters.kmers_size = 4;
    parameters.kmers_stats_fpath = "@kmers_stats_fpath";
    parameters.kmers_fpath = "@kmers_fpath";

    Pattern kmer = {1, 2, 3, 4};
    KmerIndex kmer_index = {
            {{1, 2, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {6, 6},
                                    VariantSitePath {
                                            VariantSite {5, 1}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {7, 7},
                                    VariantSitePath {
                                            VariantSite {5, 2}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {8, 8},
                                    VariantSitePath {
                                            VariantSite {5, 2}
                                    },
                                    SearchVariantSiteState::outside_variant_site
                            }
                    }
            },
            {{2, 4, 3, 4},
                    SearchStates {
                            SearchState {
                                    SA_Interval {9, 10},
                                    VariantSitePath {},
                                    SearchVariantSiteState::outside_variant_site
                            },
                            SearchState {
                                    SA_Interval {11, 11},
                                    VariantSitePath {
                                            VariantSite {5, 2},
                                            VariantSite {7, 2}
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

    sdsl::int_vector<> expected_1 = {2, 0, 2, 3, 1, 1, 1};
    sdsl::util::bit_compress(expected_1);

    sdsl::int_vector<> expected_2 = {3, 1, 1, 1, 2, 0, 2};
    sdsl::util::bit_compress(expected_2);

    bool result = stats_kmer_entry == expected_1 or stats_kmer_entry == expected_2;
    EXPECT_TRUE(result);
}
