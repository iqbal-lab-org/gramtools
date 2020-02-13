#ifndef JSON_SPEC
#define JSON_SPEC

#include <nlohmann/json.hpp>

using JSON = nlohmann::json;

namespace gram::json::spec {
    const JSON site_fields{
            {"ALS",
                    {
                            {"Desc", "Alleles at this site"}
                    }},
            {"GT",
                    {
                            {"Desc", "Sample Genotype"}
                    }},
            {"HAPG",
                    {
                            {"Desc", "Sample haplogroups of genotyped alleles"}
                    }},
            {"COVS",
                    {
                            {"Desc", "Coverage on each allele"}
                    }},
            {"DP",
                    {
                            {"Desc", "Total depth on this site"}
                    }
            }
    };

    const JSON json_prg{
            {"Model", "UNKNOWN"},
            {"Site_Fields", site_fields},
            {"Samples", JSON::object()},
            {"Sites", JSON::array()},
            {"Lvl1_Sites", JSON::array()},
            {"Child_map", JSON::object()}
    };
}

namespace gram::json{
    class Json_Prg;
    class Json_Site;
    using json_prg_ptr = std::shared_ptr<Json_Prg>;
    using json_site_ptr = std::shared_ptr<Json_Site>;
    using json_site_vec = std::vector<json_site_ptr>;

    class Json_Prg{
    protected:
        JSON json_prg;
        json_site_vec sites;
        void combine_with(const Json_Prg& other);
    public:
        Json_Prg() : json_prg(gram::json::spec::json_prg) {}
        void add_site(json_site_ptr json_site);
        JSON const& get_prg() const {return json_prg;}
        JSON get_prg_copy() const {return json_prg;}
        void set_prg(JSON const& json_prg) {this->json_prg = json_prg;}
    };

    class LevelGenotyper_Json : public Json_Prg{
    public:
        LevelGenotyper_Json() : Json_Prg() {
            json_prg.at("Model") = "LevelGenotyper";
            json_prg.at("Site_Fields")["GT_CONF"] = {
                    {"Desc", "Genotype confidence as "
                             "likelihood ratio of called and next most likely genotype."}
            };
        }
    };

    class Json_Site {
    protected:
        JSON json_site;

        void combine_with(const Json_Site &other);
        // virtual void add_model_specific_part(const Json_Site &other) = 0;

    public:
        Json_Site() {
            for (const auto &element : gram::json::spec::site_fields.items()) {
                json_site.emplace(element.key(), JSON::array());
            }
        }
        JSON get_site_copy() const {return json_site;}
        void set_site(JSON const& json_site) {this->json_site = json_site;}
    };

    class LevelGenotyped_Json_Site : public Json_Site {
    public:
        LevelGenotyped_Json_Site() : Json_Site(){
           json_site.emplace("GT_CONF", JSON::array());
        }
    };
}

#endif //JSON_SPEC
