#include "gtest/gtest.h"

#include "genotype/quasimap/coverage/grouped_allele_counts.hpp"
#include "genotype/quasimap/coverage/coverage_common.hpp"

#include "submod_resources.hpp"


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


TEST(GroupedAlleleCount, GivenSitesGroupedAlleleCounts_CorrectHashing) {
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

    // Test allele IDs in the gped allele counts are all registered and hashed
    HashSet<AlleleIds> allele_ids;
    for (auto const& entry: result) allele_ids.insert(entry.first);
    HashSet<AlleleIds> expected_allele_ids = {
            AlleleIds {1, 3},
            AlleleIds {2},
            AlleleIds {1, 4},
    };
    EXPECT_EQ(allele_ids, expected_allele_ids);

    // Test group IDs are distinct and 'full': allocated from 0 & increasing by one
    std::vector<uint64_t> group_ids;
    for (auto const& entry: result) group_ids.push_back(entry.second);
    std::sort(group_ids.begin(), group_ids.end());

    std::vector<uint64_t> expected_group_ids{0, 1, 2};
    EXPECT_EQ(group_ids, expected_group_ids);
}

TEST(GroupedAlleleCount_IDToCount, OneSite_CorrectGroupIDToCounts){
   SitesGroupedAlleleCounts sites = {
           {{AlleleIds {0, 1}, 19}, {AlleleIds {0}, 2}}
   };

    AlleleGroupHash allele_ids_groups_hash = {
            {AlleleIds{0},    0},
            {AlleleIds{0, 1}, 1},
    };

    SitesGroupIDToCounts expected{ // Ordered by the key
            {{"0", 2}, {"1", 19}}
    };
    auto result = get_group_id_counts(sites, allele_ids_groups_hash);
    EXPECT_EQ(result, expected);
}

TEST(GroupedAlleleCount_IDToCount, TwoSites_CorrectGroupIDToCounts){
    SitesGroupedAlleleCounts sites = {
            {
                    {AlleleIds {1, 3}, 1},
                    {AlleleIds {1, 4}, 2}
            },
            {
                    {AlleleIds {2}, 10},
                    {AlleleIds {3, 4}, 2},
                    {AlleleIds {1, 3}, 20},
            },
    };
    AlleleGroupHash allele_ids_groups_hash = {
            {AlleleIds {1, 3}, 0},
            {AlleleIds {1, 4}, 1},
            {AlleleIds {2}, 2},
            {AlleleIds {3, 4}, 3},
    };

    SitesGroupIDToCounts expected{
            {
                    {"0", 1},
                    {"1", 2}
            },
            {
                    {"0", 20},
                    {"2", 10},
                    {"3", 2},
            },
    };
    auto result = get_group_id_counts(sites, allele_ids_groups_hash);
    EXPECT_EQ(result, expected);
}

TEST(ReverseAlleleGroupHash, Succeeds){
    AlleleGroupHash allele_ids_groups_hash = {
            {AlleleIds {1, 3}, 42},
            {AlleleIds {1, 4}, 43}
    };
    GroupIDToAlleles reversed = get_group_id_alleles(allele_ids_groups_hash);
    GroupIDToAlleles expected{
            {"42", AlleleIds{1, 3}},
            {"43", AlleleIds{1, 4}},
    };
    EXPECT_EQ(reversed, expected);
}

TEST(ReverseAlleleGroupHash, IsOrderedByNumericValue){
    // Lexicographically, "30" < "9", but gets ordered by numeric value
    AlleleGroupHash allele_ids_groups_hash = {
            {AlleleIds {1, 3}, 30},
            {AlleleIds {1, 4}, 9}
    };
    GroupIDToAlleles reversed = get_group_id_alleles(allele_ids_groups_hash);
    GroupIDToAlleles expected{
            {"9", AlleleIds{1, 4}},
            {"30", AlleleIds{1, 3}},
    };
    EXPECT_EQ(reversed, expected);
}

/*
 * The underlying data structures are maps, such that one can always expect
 * the groups to be listed in increasing integer order, both for `site_counts` (for each site)
 * and for `allele_groups`.
 */
class TestGetJSON : public ::testing::Test{
protected:
    void SetUp(){
        site1 = {
                {AlleleIds {1, 3}, 1},
                {AlleleIds {1, 4}, 2}
        };
        site2 = {
                {AlleleIds {0}, 19},
                {AlleleIds {1, 4}, 5}
        };

        group_ids = AlleleGroupHash{
                {AlleleIds{1, 3}, 0},
                {AlleleIds{1, 4}, 2},
                {AlleleIds{0}, 1},
        };
    }
    SitesGroupedAlleleCounts sites;
    GroupedAlleleCounts site1, site2;
    AlleleGroupHash group_ids;
    std::string expected_allele_groups = R"({"0":[1,3],"1":[0],"2":[1,4]})";
    std::string expected_site_one_counts = R"([{"0":1,"2":2}])";
    std::string expected_site_two_counts = R"([{"1":19,"2":5}])";
    std::string expected_all_counts = R"([{"0":1,"2":2},{"1":19,"2":5}])";
};

TEST_F(TestGetJSON, AlleleIds_CorrectJson) {
    auto result = get_json(sites, group_ids);
    auto result_allele_groups =
            result.at("grouped_allele_counts").at("allele_groups").dump();
    EXPECT_EQ(result_allele_groups, expected_allele_groups);
}


TEST_F(TestGetJSON, SiteOne_CorrectJson) {
    sites.push_back(site1);
    auto result = get_json(sites, group_ids);
    auto result_site_counts = result.at("grouped_allele_counts").at("site_counts").dump();
    EXPECT_EQ(result_site_counts, expected_site_one_counts);
}

TEST_F(TestGetJSON, SiteTwo_CorrectJson) {
    sites.push_back(site2);
    auto result = get_json(sites, group_ids);
    auto result_site_counts = result.at("grouped_allele_counts").at("site_counts").dump();
    EXPECT_EQ(result_site_counts, expected_site_two_counts);
}

// Empty (no coverage) sites get an empty entry
TEST_F(TestGetJSON, EmptySites_CorrectJson) {
    sites.push_back(GroupedAlleleCounts{});
    sites.push_back(GroupedAlleleCounts{});
    auto result = get_json(sites, group_ids);
    auto result_site_counts = result.at("grouped_allele_counts").at("site_counts").dump();
    std::string expected_site_counts = R"([{},{}])";
    EXPECT_EQ(result_site_counts, expected_site_counts);
}


TEST_F(TestGetJSON, TwoSites_CorrectFullJson) {
    sites.push_back(site1);
    sites.push_back(site2);
    auto result = get_json(sites, group_ids).dump();
    // Entries get alphabetically sorted
    std::string expected = std::string(R"({"grouped_allele_counts":{"allele_groups":)")
            + expected_allele_groups
            + std::string(R"(,"site_counts":)")
            + expected_all_counts
            + std::string("}}");
    EXPECT_EQ(result, expected);
}

