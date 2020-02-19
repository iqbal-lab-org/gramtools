#ifndef PRG_JSON_SPEC
#define PRG_JSON_SPEC

#include "common.hpp"

namespace gram::json{

    class Json_Prg{
    protected:
        JSON json_prg;
        json_site_vec sites;
    public:
        Json_Prg() : json_prg(gram::json::spec::json_prg) {}
        void add_samples(const Json_Prg &other, const bool force = false);
        void combine_with(const Json_Prg &other, bool force = false);
        void set_sample_info(std::string const& name, std::string const& desc);

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

    class Simulated_Json : public Json_Prg{
    public:
        Simulated_Json() : Json_Prg() {
            json_prg.at("Model") = "Simulated_Path";
        }
    };
}

#endif //PRG_JSON_SPEC