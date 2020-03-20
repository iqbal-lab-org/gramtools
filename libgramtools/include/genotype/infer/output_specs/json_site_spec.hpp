#ifndef SITE_JSON_SPEC
#define SITE_JSON_SPEC

#include "common/utils.hpp"
#include "json_common.hpp"

namespace gram::json {
    struct site_rescaler {
        std::size_t index;
        AlleleId hapg;
    };
    using allele_combi_map = std::map<std::string, site_rescaler>;
    using allele_vec = std::vector<std::string>;

    class Json_Site {
    protected:
        JSON json_site;

    public:
        Json_Site() {
            auto site_fields = spec::json_site_fields();
            for (const auto &element : site_fields.items()) {
                json_site.emplace(element.key(), JSON::array());
            }
            json_site.emplace("SEG", "");
        }

        Json_Site(Json_Site const &other) : json_site(other.json_site) {}

        // Functions implementing site combining
        void build_allele_combi_map(JSON const &json_site, allele_combi_map &m);

        void append_entries_from(JSON const &json_site);

        allele_vec get_all_alleles(allele_combi_map &m);

        JSON rescale_entries(allele_combi_map const &m) const;

        void combine_with(const Json_Site &other);

        std::size_t get_pos() const { return json_site.at("POS"); }
        JSON const &get_site() const { return json_site; }
        JSON get_site_copy() const { return json_site; }

        void set_site(JSON const &json_site) { this->json_site = json_site; }
        void set_segment(std::string ID){ this->json_site.at("SEG") = ID; }
    };

    class LevelGenotyped_Json_Site : public Json_Site {
    public:
        LevelGenotyped_Json_Site() : Json_Site() {
            json_site.emplace("GT_CONF", JSON::array());
        }
    };
}

#endif //SITE_JSON_SPEC
