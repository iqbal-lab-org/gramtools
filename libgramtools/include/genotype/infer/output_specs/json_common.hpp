#ifndef COMMON_JSON_SPEC
#define COMMON_JSON_SPEC

#include <nlohmann/json.hpp>
#include "fields.hpp"

using JSON = nlohmann::json;
using namespace gram::genotype::output_spec;

namespace gram::json{
    class Json_Prg;
    class Json_Site;
    using json_prg_ptr = std::shared_ptr<Json_Prg>;
    using json_site_ptr = std::shared_ptr<Json_Site>;
    using json_site_vec = std::vector<json_site_ptr>;

    class JSONParseException : public std::exception{
    protected:
        std::string msg;
    public:
        virtual char const* what() const throw(){ return msg.c_str(); }
    };

    class JSONCombineException : public JSONParseException{
    public:
        JSONCombineException(std::string msg) {
            this->msg = std::string("JSONCombineException: "+ msg); }
    };

    class JSONConsistencyException : public JSONParseException{
    public:
        JSONConsistencyException(std::string msg) {
            this->msg = std::string("JSONConsistencyException: "+ msg); }
    };
}

namespace gram::json::spec {

    static JSON json_site_fields(){
        header_vec h = vcf_format_headers();
        JSON result =
                {
                        {"ALS",
                                {{"Desc", "Alleles at this site"}}
                        },
                        {"HAPG",
                                {{"Desc", "Sample haplogroups of genotyped alleles"}},
                        }
                };

        // Populate with headers common with vcf output
        for (auto const& header : h){
            result[header.ID] = { {"Desc", header.desc} };
        }
        return result;
    }

    const JSON json_prg{
            {"Model", "UNKNOWN"},
            {"Site_Fields", json_site_fields()},
            {"Samples", JSON::array()},
            {"Sites", JSON::array()},
            {"Lvl1_Sites", JSON::array()},
            {"Child_Map", JSON::object()}
    };
}

#endif //COMMON_JSON_SPEC
