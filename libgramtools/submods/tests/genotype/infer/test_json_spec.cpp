#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "genotype/infer/json_spec.hpp"

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
    MockJsonSite(const MockJsonSite& other){
        this->json_site = other.json_site;
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

    MOCK_METHOD(void, add_model_specific_part, (Json_Site const&), (override));
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

public:
    json_site_vec site1_samples;
   JSON_data_store(){
        set_site1_samples();
   }
};

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

TEST_F(PRG_Combine_Fail, GivenSameJSONs_DoesNotFail){
    json_prg2.set_prg(the_prg);
    ASSERT_NO_THROW(json_prg1.combine_with(json_prg2));
}

TEST_F(PRG_Combine_Fail, GivenDifferentModels_Fails){
    the_prg.at("Model") = "A_different_model";
    json_prg2.set_prg(the_prg);
    ASSERT_THROW(json_prg1.combine_with(json_prg2), JSONCombineException);
}

TEST_F(PRG_Combine_Fail, GivenDifferentPRGs_Fails){
    auto copy = the_prg;
    the_prg.at("Child_Map") = {};
    json_prg2.set_prg(the_prg);
    EXPECT_THROW(json_prg1.combine_with(json_prg2), JSONCombineException);

    EXPECT_EQ(copy, json_prg1.get_prg());
    copy.at("Lvl1_Sites").push_back("all");
    json_prg2.set_prg(copy);
    EXPECT_THROW(json_prg1.combine_with(json_prg2), JSONCombineException);
}

TEST_F(PRG_Combine_Fail, GivenDifferentSiteSpecs_Fails){
    the_prg.at("Site_Fields").at("GT").at("Desc") = "Greater Than";
    json_prg2.set_prg(the_prg);
    ASSERT_THROW(json_prg1.combine_with(json_prg2), JSONCombineException);
}

TEST_F(PRG_Combine_Fail, GivenDifferentNumOfSites_Fails){
    json_prg2.set_prg(the_prg);
    json_prg2.add_site(std::make_shared<LevelGenotyped_Json_Site>());
    ASSERT_THROW(json_prg1.combine_with(json_prg2), JSONCombineException);
}

class Site_Combine_Fail : public ::testing::Test{
protected:
   void SetUp(){
       JSON_data_store data;
       the_site_json = data.site1_samples.at(0)->get_site_copy();
       // It this site: MockJsonSite sample1({"CTCCT", "CTT"}, {0, 0}, {0, 0}, {10, 2}, {11});
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

TEST(Site_Json_AppendEntries, GivenTwoSites_CorrectAppending){
    //MockJsonSite sample1({"CTCCT", "CTT"}, {0, 0}, {0, 0}, {10, 2}, {11});
    //MockJsonSite sample2({"CTCCT", "CTT"}, {1, 1}, {1, 1}, {2, 10}, {11});
    JSON_data_store data;
    auto sample1 = data.site1_samples.at(0), sample2 = data.site1_samples.at(1);
    sample1->append_entries_from(sample2->get_site());

    MockJsonSite expected({"CTCCT", "CTT"}, {{0, 0}, {1, 1}},
            {{0, 0}, {1, 1}}, {{10, 2}, {2, 10}}, {11, 11});
    EXPECT_EQ(sample1->get_site(), expected.get_site());
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
