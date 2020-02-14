#ifndef JSON_SPEC
#define JSON_SPEC

#include <nlohmann/json.hpp>
#include <iostream>
#include "common/utils.hpp"

using JSON = nlohmann::json;

namespace gram::genotype::infer{
    using GtypedIndex = std::size_t; /**< The index of an allele in an allele vector */
    using GtypedIndices = std::vector<GtypedIndex>;
    using allele_coverages = std::vector<double>;
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
            {"Samples", JSON::object()},
            {"Sites", JSON::array()},
            {"Lvl1_Sites", JSON::array()},
            {"Child_Map", JSON::object()}
    };
}

namespace gram::json{

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

    class Json_Prg;
    class Json_Site;
    using json_prg_ptr = std::shared_ptr<Json_Prg>;
    using json_site_ptr = std::shared_ptr<Json_Site>;
    using json_site_vec = std::vector<json_site_ptr>;

    class Json_Prg{
    protected:
        JSON json_prg;
        json_site_vec sites;
    public:
        Json_Prg() : json_prg(gram::json::spec::json_prg) {}
        void combine_with(const Json_Prg& other);

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

    struct site_rescaler{
        std::size_t index;
        AlleleId hapg;
    };
    using allele_combi_map = std::map<std::string, site_rescaler>;
    using allele_vec = std::vector<std::string>;

    class Json_Site {
    protected:
        JSON json_site;

        virtual void add_model_specific_part(const Json_Site &other) = 0;

    public:
        Json_Site() {
            for (const auto &element : gram::json::spec::site_fields.items()) {
                json_site.emplace(element.key(), JSON::array());
            }
        }
        Json_Site(Json_Site const& other) : json_site(other.json_site) {}

        // Functions implementing site combining
        void build_allele_combi_map(JSON const& json_site, allele_combi_map& m);
        void append_entries_from(JSON const& json_site);
        allele_vec get_all_alleles(allele_combi_map& m);
        JSON rescale_entries(allele_combi_map const &m) const;
        void combine_with(const Json_Site &other);

        JSON const& get_site() const {return json_site;}
        JSON get_site_copy() const {return json_site;}
        void set_site(JSON const& json_site) {this->json_site = json_site;}
    };

    class LevelGenotyped_Json_Site : public Json_Site {
    public:
        LevelGenotyped_Json_Site() : Json_Site(){
           json_site.emplace("GT_CONF", JSON::array());
        }

        void add_model_specific_part(const Json_Site &other) override {};
    };
}

#endif //JSON_SPEC
