#include "gtest/gtest.h"

#include "genotype/quasimap/coverage/grouped_allele_counts.hpp"
#include "genotype/quasimap/coverage/common.hpp"
#include "src_common/common.hpp"


TEST(GroupedAlleleCount, GivenTwoVariantSites_CorrectEmptySitesVectorSize) {
    auto prg_raw = encode_prg("gct5c6g6t6ac7cc8a8");
    auto prg_info = generate_prg_info(prg_raw);
    auto grouped_allele_counts = coverage::generate::grouped_allele_counts(prg_info);

    auto result = grouped_allele_counts.size();
    uint64_t expected = 2;
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, GivenTwoSearchStates_CorrectCoverage) {
    auto prg_raw = encode_prg("gct5c6g6t6ac7cc8a8");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uniqueLoci compatible_loci = {
            VariantLocus{7, FIRST_ALLELE},
            VariantLocus{5, FIRST_ALLELE},
            VariantLocus{5, FIRST_ALLELE + 1}

    };
    coverage::record::grouped_allele_counts(coverage, compatible_loci);
    auto result = coverage.grouped_allele_counts;
    SitesGroupedAlleleCounts expected = {
            GroupedAlleleCounts {{AlleleIds {0, 1}, 1}},
            GroupedAlleleCounts {{AlleleIds {0}, 1}},
    };
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, GivenSingleSearchState_CorrectCoverage) {
    auto prg_raw = encode_prg("gct5c6g6t6ac7cc8a8");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uniqueLoci compatible_loci = {
            VariantLocus{5, FIRST_ALLELE + 2},

    };
    coverage::record::grouped_allele_counts(coverage, compatible_loci);
    auto result = coverage.grouped_allele_counts;
    SitesGroupedAlleleCounts expected = {
            GroupedAlleleCounts {{AlleleIds {2}, 1}},
            GroupedAlleleCounts {},
    };
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, MultipleReads_CorrectCoverage) {
    auto prg_raw = encode_prg("gct5c6g6t6ac7cc8a8");
    auto prg_info = generate_prg_info(prg_raw);
    auto coverage = coverage::generate::empty_structure(prg_info);

    uniqueLoci read1_compatible_loci = {
            VariantLocus{7, FIRST_ALLELE + 1},
            VariantLocus{5, FIRST_ALLELE + 2},
            VariantLocus{5, FIRST_ALLELE}

    };
    uniqueLoci read2_compatible_loci = {
            VariantLocus{7, FIRST_ALLELE + 1},
            VariantLocus{5, FIRST_ALLELE + 3},
            VariantLocus{5, FIRST_ALLELE}

    };

    coverage::record::grouped_allele_counts(coverage,
                                            read1_compatible_loci);
    coverage::record::grouped_allele_counts(coverage,
                                            read2_compatible_loci);

    auto result = coverage.grouped_allele_counts;
    SitesGroupedAlleleCounts expected = {
            GroupedAlleleCounts {
                    {AlleleIds {0, 2}, 1},
                    {AlleleIds {0, 3}, 1}
            },
            GroupedAlleleCounts {
                    {AlleleIds {1}, 2}
            }
    };
    EXPECT_EQ(result, expected);
}


bool valid_hash_allele_groups(const AlleleGroupHash &allele_ids_groups_hash,
                              const HashSet<AlleleIds> &correct_allele_ids_groups) {
    std::unordered_set<uint64_t> seen_hashes;
    for (const auto &entry: allele_ids_groups_hash) {
        auto allele_ids = entry.first;
        bool allele_ids_correct = correct_allele_ids_groups.find(allele_ids)
                                  != correct_allele_ids_groups.end();
        if (not allele_ids_correct)
            return false;

        auto hash = entry.second;
        bool hash_seen = seen_hashes.find(hash) != seen_hashes.end();
        if (hash_seen)
            return false;
        seen_hashes.insert(hash);
    }
    return true;
}


TEST(GroupedAlleleCount, GivenSitesGroupedAlleleCounts_CorrectlyAssignHashValuesToAlleleIdsGroups) {
    SitesGroupedAlleleCounts grouped_allele_counts = {
            GroupedAlleleCounts {
                    {AlleleIds {1, 3}, 1},
                    {AlleleIds {1, 4}, 1}
            },
            GroupedAlleleCounts {
                    {AlleleIds {2}, 2}
            }
    };
    auto result = hash_allele_groups(grouped_allele_counts);
    HashSet<AlleleIds> expected = {
            AlleleIds {1, 3},
            AlleleIds {1, 4},
            AlleleIds {2},
    };
    auto correct = valid_hash_allele_groups(result, expected);
    EXPECT_TRUE(correct);
}



TEST(GroupedAlleleCount, GivenSingleSite_CorrectJsonString) {
    GroupedAlleleCounts site = {
            {AlleleIds {1, 3}, 1},
            {AlleleIds {1, 4}, 2}
    };
    AlleleGroupHash allele_ids_groups_hash = {
            {AlleleIds {1, 3}, 42},
            {AlleleIds {1, 4}, 43}
    };
    auto result = dump_site(allele_ids_groups_hash, site);
    std::string expected = R"({"43":2,"42":1})";
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, GivenMultipleSites_CorrectSitesCountsJsonString) {
    SitesGroupedAlleleCounts sites = {
            GroupedAlleleCounts {
                    {AlleleIds {1, 3}, 1},
                    {AlleleIds {1, 4}, 3}
            },
            GroupedAlleleCounts {
                    {AlleleIds {2}, 2}
            }
    };
    AlleleGroupHash allele_ids_groups_hash = {
            {AlleleIds {1, 3}, 42},
            {AlleleIds {1, 4}, 43},
            {AlleleIds {2},    44}
    };
    auto result = dump_site_counts(allele_ids_groups_hash, sites);
    std::string expected = R"("site_counts":[{"43":3,"42":1},{"44":2}])";
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, GivenHashedAlleleIdsGroups_CorrectAlleleGroupsJsonString) {
    AlleleGroupHash allele_ids_groups_hash = {
            {AlleleIds {1, 3}, 42},
            {AlleleIds {1, 4}, 43},
            {AlleleIds {2},    44}
    };
    auto result = dump_allele_groups(allele_ids_groups_hash);
    std::string expected = R"("allele_groups":{"44":[2],"43":[1,4],"42":[1,3]})";
    EXPECT_EQ(result, expected);
}


TEST(GroupedAlleleCount, GivenMultipleSites_CorrectFullJsonString) {
    SitesGroupedAlleleCounts sites = {
            GroupedAlleleCounts {
                    {AlleleIds {1, 3}, 1},
                    {AlleleIds {1, 4}, 3}
            },
            GroupedAlleleCounts {
                    {AlleleIds {2}, 2}
            }
    };
    auto result = dump_grouped_allele_counts(sites);
    std::string expected = R"({"grouped_allele_counts":{"site_counts":[{"0":3,"1":1},{"2":2}],"allele_groups":{"0":[1,4],"2":[2],"1":[1,3]}}})";
    std::cout << expected << std::endl;
    EXPECT_EQ(result, expected);
}
