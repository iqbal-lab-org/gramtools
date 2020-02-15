#ifndef COMMON_JSON_SPEC
#define COMMON_JSON_SPEC

#include <nlohmann/json.hpp>

using JSON = nlohmann::json;

namespace gram::genotype::infer{
    using GtypedIndex = std::size_t; /**< The index of an allele in an allele vector */
    using GtypedIndices = std::vector<GtypedIndex>;
    using allele_coverages = std::vector<double>;
}

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
            {"Samples", JSON::array()},
            {"Sites", JSON::array()},
            {"Lvl1_Sites", JSON::array()},
            {"Child_Map", JSON::object()}
    };
}

#endif //COMMON_JSON_SPEC
