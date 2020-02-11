#ifndef JSON_SPEC
#define JSON_SPEC

#include <nlohmann/json.hpp>

using JSON = nlohmann::json;

namespace json_::spec {
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

#endif //JSON_SPEC
