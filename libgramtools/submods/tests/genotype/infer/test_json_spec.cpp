#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "genotype/infer/json_spec.hpp"

using namespace gram::json;

class MockJsonSite : public Json_Site {
public:
    MOCK_METHOD(void, add_model_specific_part, (Json_Site const&), (override));
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
       the_site_json = fixed_site.get_site_copy();
       the_site_json.at("ALS") = JSON::array({"CTCCT", "CTT"});
       the_site_json.at("GT").push_back(JSON::array({0, 0}));
       the_site_json.at("HAPG").push_back(JSON::array({0, 0}));
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
