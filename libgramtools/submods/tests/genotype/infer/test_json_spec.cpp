#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "genotype/infer/json_spec/prg_spec.hpp"
#include "genotype/infer/json_spec/site_spec.hpp"

using namespace gram;
using namespace gram::json;
using namespace gram::genotype::infer;

namespace gram::json{
    bool operator==(site_rescaler const& first, site_rescaler const& second){
        return (first.index == second.index) && (first.hapg == second.hapg);
    }
}

class MockJsonSite : public Json_Site {
public:
    MockJsonSite() : Json_Site() {}
    MockJsonSite(MockJsonSite const& other){
        this->json_site = other.json_site;
    }

    void set_json(JSON const& json){ this->json_site = json;}

    void make_null(){
        if (json_site.at("GT").size() != 1) return;
        json_site.at("GT").at(0) = JSON::array({nullptr});
    }

    MockJsonSite(allele_vec als, GtypedIndices gts, AlleleIds hapgs,
                 allele_coverages coverages, std::size_t total_cov) : Json_Site(){
        json_site.at("ALS") = JSON(als);
        json_site.at("GT").push_back(JSON(gts));
        json_site.at("HAPG").push_back(JSON(hapgs));
        json_site.at("COVS").push_back(JSON(coverages));
        json_site.at("DP").push_back(total_cov);
    }

    MockJsonSite(allele_vec als, std::vector<GtypedIndices> gts, std::vector<AlleleIds> hapgs,
                 std::vector<allele_coverages> coverages, std::vector<std::size_t> total_covs) : Json_Site(){
        json_site.at("ALS") = JSON(als);
        for(auto entry : gts) json_site.at("GT").push_back(JSON(entry));
        for(auto entry : hapgs) json_site.at("HAPG").push_back(JSON(entry));
        for(auto entry : coverages) json_site.at("COVS").push_back(JSON(entry));
        for(auto entry : total_covs) json_site.at("DP").push_back(JSON(entry));
    }
};

class JSON_data_store{
private:
    void set_site1_samples(){
        MockJsonSite sample1({"CTCCT", "CTT"}, {0, 0}, {0, 0}, {10, 2}, {11});
        site1_samples.push_back(std::make_shared<MockJsonSite>(sample1));

        MockJsonSite sample2({"CTCCT", "CTT"}, {1, 1}, {1, 1}, {2, 10}, {11});
        site1_samples.push_back(std::make_shared<MockJsonSite>(sample2));

        MockJsonSite sample3({"CTCCT", "GTT"}, {0, 1}, {0, 2}, {5, 5}, {12});
        site1_samples.push_back(std::make_shared<MockJsonSite>(sample3));
    }

    void set_site2_samples() {
        MockJsonSite sample1({"AAAAAAA", "AAA"}, {1}, {1}, {20, 1}, {23});
        site2_samples.push_back(std::make_shared<MockJsonSite>(sample1));

        MockJsonSite sample2({"AAAAAAA", "A"}, {1}, {4}, {0, 18}, {24});
        site2_samples.push_back(std::make_shared<MockJsonSite>(sample2));
    }

    void set_prgs(){
        prg1.set_sample_info("Gazorp", "");
        prg1.add_site(site1_samples.at(0));
        prg1.add_site(site2_samples.at(0));

        prg2.set_sample_info("Dorp", "");
        prg2.add_site(site1_samples.at(1));
        prg2.add_site(site2_samples.at(1));
    }

public:
    json_site_vec site1_samples, site2_samples;
    Json_Prg prg1, prg2;
    JSON_data_store(){
        set_site1_samples();
        set_site2_samples();
        set_prgs();
    }
};


class Site_Combine_Fail : public ::testing::Test{
protected:
    void SetUp(){
        JSON_data_store data;
        the_site_json = data.site1_samples.at(0)->get_site_copy();
        // site: MockJsonSite sample1({"CTCCT", "CTT"}, {0, 0}, {0, 0}, {10, 2}, {11});
        fixed_site.set_site(the_site_json);
    }
    JSON the_site_json;
    MockJsonSite fixed_site, test_site;
};

TEST_F(Site_Combine_Fail, GivenSameJSONs_NoFail){
    test_site.set_site(the_site_json);
    ASSERT_NO_THROW(fixed_site.combine_with(test_site));
}

TEST_F(Site_Combine_Fail, GivenDifferentREFAllele_Fails){
    the_site_json.at("ALS").at(0) = "NOTSAME";
    test_site.set_site(the_site_json);
    ASSERT_THROW(fixed_site.combine_with(test_site), JSONCombineException);
}

TEST_F(Site_Combine_Fail, GivenInconsistenHAPGs_Fails){
    the_site_json.at("HAPG").at(0).at(0) = 1;
    test_site.set_site(the_site_json);
    ASSERT_THROW(fixed_site.combine_with(test_site), JSONConsistencyException);
}

TEST_F(Site_Combine_Fail, GivenDifferentCOVandALSCardinality_Fails){
    the_site_json.at("COVS").at(0) = JSON::array({10});
    test_site.set_site(the_site_json);
    ASSERT_THROW(fixed_site.combine_with(test_site), JSONConsistencyException);
}

TEST(Site_Json_CombiMap, Add2Samples_CorrectCombiMap){
    JSON_data_store data;
    allele_combi_map result;
    MockJsonSite site;

    JSON sample1 = data.site1_samples.at(0)->get_site_copy();
    site.build_allele_combi_map(sample1, result);

    JSON sample2 = data.site1_samples.at(1)->get_site_copy();
    site.build_allele_combi_map(sample2, result);

    allele_combi_map expected{
            {"CTCCT", site_rescaler{0, 0}},
            {"CTT", site_rescaler{1, 1}},
    };
    EXPECT_EQ(result, expected);
}

TEST(Site_Json_RescaleEntries, GivenCombiMap_CorrectRescaledJSON){
    //MockJsonSite sample2({"CTCCT", "CTT"}, {1, 1}, {1, 1}, {2, 10}, {11});
    allele_combi_map m{
            {"CTCCT", site_rescaler{0, 0}},
            {"CCC", site_rescaler{1, 2}},
            {"CTT", site_rescaler{2, 1}},
    };

    JSON_data_store data;
    auto sample2 = data.site1_samples.at(1);
    auto result = sample2->rescale_entries(m);
    MockJsonSite expected_site({"CTCCT", "CTT"}, {2, 2}, {1, 1}, {2, 0, 10}, {11});
    EXPECT_EQ(result, expected_site.get_site());
}

TEST(Site_Json_AppendEntries, GivenTwoGtedSites_CorrectAppending){
    //MockJsonSite sample1({"CTCCT", "CTT"}, {0, 0}, {0, 0}, {10, 2}, {11});
    //MockJsonSite sample2({"CTCCT", "CTT"}, {1, 1}, {1, 1}, {2, 10}, {11});
    JSON_data_store data;
    auto sample1 = data.site1_samples.at(0), sample2 = data.site1_samples.at(1);

    sample1->combine_with(*sample2);
    MockJsonSite expected({"CTCCT", "CTT"}, {{0, 0}, {1, 1}},
                          {{0, 0}, {1, 1}}, {{10, 2}, {2, 10}}, {11, 11});
    EXPECT_EQ(sample1->get_site(), expected.get_site());
}

TEST(Site_Combine_Success, GivenOneNullGtSite_Succeeds){
    //MockJsonSite sample1({"CTCCT", "CTT"}, {0, 0}, {0, 0}, {10, 2}, {11});
    JSON_data_store data;
    auto sample1 = data.site1_samples.at(0);
    MockJsonSite to_null_site;
    to_null_site.set_json(sample1->get_site());
    to_null_site.make_null();
    auto sample2 = std::make_shared<MockJsonSite>(to_null_site);
    sample1->combine_with(*sample2);

    auto json_result = sample1->get_site();
    auto expected_gt_first = GtypedIndices{0, 0};
    EXPECT_EQ(json_result.at("GT").at(0), expected_gt_first);
    EXPECT_EQ(json_result.at("GT").at(1), JSON::array({nullptr}));
}


TEST(Site_Combine_Success, GivenThreeSites_CorrectCombinedSite){
    //MockJsonSite sample1({"CTCCT", "CTT"}, {0, 0}, {0, 0}, {10, 2}, {11});
    //MockJsonSite sample2({"CTCCT", "CTT"}, {1, 1}, {1, 1}, {2, 10}, {11});
    //MockJsonSite sample3({"CTCCT", "GTT"}, {0, 1}, {0, 2}, {5, 5}, {12});
    JSON_data_store data;
    auto sample1 = data.site1_samples.at(0), sample2 = data.site1_samples.at(1),
            sample3 = data.site1_samples.at(2);
    sample1->combine_with(*sample2);
    sample1->combine_with(*sample3);

    MockJsonSite expected({"CTCCT", "CTT", "GTT"},
                          {{0, 0}, {1, 1}, {0, 2}},
                          {{0, 0}, {1, 1}, {0, 2}},
                          {{10, 2, 0}, {2, 10, 0}, {5, 0, 5}},
                          {11, 11, 12});
    EXPECT_EQ(sample1->get_site(), expected.get_site());

    // Now show associativity: (sample1 + sample2) + sample 3 == sample1 + (sample2 + sample3)
    JSON_data_store data_again;
    sample1 = data_again.site1_samples.at(0), sample2 = data_again.site1_samples.at(1),
    sample3 = data_again.site1_samples.at(2);
    sample2->combine_with(*sample3);
    sample1->combine_with(*sample2);
    EXPECT_EQ(sample1->get_site(), expected.get_site());
}

class PRG_Combine_Fail : public ::testing::Test{
protected:
    void SetUp(){
        the_prg = gram::json::spec::json_prg;
        the_prg.at("Model") = "M1";
        the_prg.at("Child_Map") = {
                {0,{
                           {1,JSON::array({2,3})}
                   }}
        };
        the_prg.at("Lvl1_Sites").push_back(0);

        json_prg1.set_prg(the_prg);
    }
    JSON the_prg;
    Json_Prg json_prg1;
    Json_Prg json_prg2;
};


TEST_F(PRG_Combine_Fail, GivenDifferentModels_Fails){
    the_prg.at("Model") = "A_different_model";
    json_prg2.set_prg(the_prg);
    ASSERT_THROW(json_prg1.combine_with(json_prg2, false), JSONCombineException);
}

TEST_F(PRG_Combine_Fail, GivenDifferentPRGs_Fails){
    auto copy = the_prg;
    the_prg.at("Child_Map") = {};
    json_prg2.set_prg(the_prg);
    EXPECT_THROW(json_prg1.combine_with(json_prg2, false), JSONCombineException);

    EXPECT_EQ(copy, json_prg1.get_prg());
    copy.at("Lvl1_Sites").push_back("all");
    json_prg2.set_prg(copy);
    EXPECT_THROW(json_prg1.combine_with(json_prg2, false), JSONCombineException);
}

TEST_F(PRG_Combine_Fail, GivenDifferentSiteSpecs_Fails){
    the_prg.at("Site_Fields").at("GT").at("Desc") = "Greater Than";
    json_prg2.set_prg(the_prg);
    ASSERT_THROW(json_prg1.combine_with(json_prg2, false), JSONCombineException);
}

TEST_F(PRG_Combine_Fail, GivenDifferentNumOfSites_Fails){
    json_prg2.set_prg(the_prg);
    json_prg2.add_site(std::make_shared<LevelGenotyped_Json_Site>());
    ASSERT_THROW(json_prg1.combine_with(json_prg2, false), JSONCombineException);
}

TEST(PRG_Combine_SampleNames, GivenNamedJSONs_CanForceOrNotForceMerge){
    JSON_data_store data;
    Json_Prg prg1, prg2;
    prg1.add_site(data.site1_samples.at(0));
    prg2.add_site(data.site1_samples.at(1));
    prg1.set_sample_info("Sample1", "I am sample1");
    prg2.set_sample_info("Sample1", "I am another sample but I was named the same. Sorry.");

    EXPECT_THROW(prg1.add_samples(prg2), JSONConsistencyException);

    JSON expected = JSON::array();
    expected.push_back(prg1.get_prg().at("Samples").at(0));
    expected.push_back(prg2.get_prg().at("Samples").at(0));
    expected.at(1).at("Name") = "Sample1_1";
    prg1.add_samples(prg2, true); // Forcing duplicate samp names to be allowed
    EXPECT_EQ(prg1.get_prg().at("Samples"), expected);
}

TEST(PRG_Combine_Success, GivenTwoPrgs_CorrectCombined){
    JSON_data_store data;

    MockJsonSite s1_cpy, s2_cpy;
    s1_cpy.set_site(data.site1_samples.at(0)->get_site());
    s2_cpy.set_site(data.site2_samples.at(0)->get_site());
    data.prg1.combine_with(data.prg2);

    s1_cpy.combine_with(*data.site1_samples.at(1));
    EXPECT_EQ(data.prg1.get_prg().at("Sites").at(0), s1_cpy.get_site());

    s2_cpy.combine_with(*data.site2_samples.at(1));
    EXPECT_EQ(data.prg1.get_prg().at("Sites").at(1), s2_cpy.get_site());
}
