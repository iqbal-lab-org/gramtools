#ifndef COMMON_JSON_SPEC
#define COMMON_JSON_SPEC

#include <nlohmann/json.hpp>
#include "fields.hpp"

using JSON = nlohmann::json;
using namespace gram::genotype::output_spec;

namespace gram::json {
    class Json_Prg;

    class Json_Site;

    using json_prg_ptr = std::shared_ptr<Json_Prg>;
    using json_site_ptr = std::shared_ptr<Json_Site>;
    using json_site_vec = std::vector<json_site_ptr>;


    class JSONCombineException : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    class JSONConsistencyException : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };
}

namespace gram::json::spec {

    static JSON site_fields(){
        JSON result =
                {
                        {"POS",
                                {{"Desc", "Position on reference or pseudo-reference"}}
                        },
                        {"SEG",
                                {{"Desc", "Segment ID"}}
                        },
                        {"ALS",
                                {{"Desc", "Alleles at this site"}}
                        },
                        {"HAPG",
                                {{"Desc", "Sample haplogroups of genotyped alleles"}},
                        },
                };

        // Populate with headers common with vcf output
        for (auto const& header : common_headers()){
            if (header.meta_type == "FORMAT")
                result[header.ID] = { {"Desc", header.desc} };
        }
        return result;
    }

    static JSON filters(){
        JSON result;
        // Populate with headers common with vcf output
        for (auto const& header : common_headers()){
            if (header.meta_type == "FILTER")
                result[header.ID] = { {"Desc", header.desc} };
        }
        return result;
    }

    const JSON json_prg{
            {"Model",       "UNKNOWN"},
            {"Site_Fields", site_fields()},
            {"Filters",     filters()},
            {"Samples",     JSON::array()},
            {"Sites",       JSON::array()},
            {"Lvl1_Sites",  JSON::array()},
            {"Child_Map",   JSON::object()}
    };
}

#endif //COMMON_JSON_SPEC
