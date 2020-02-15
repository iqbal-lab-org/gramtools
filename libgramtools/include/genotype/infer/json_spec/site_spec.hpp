#ifndef SITE_JSON_SPEC
#define SITE_JSON_SPEC

#include "common.hpp"
#include "common/utils.hpp"

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

        virtual void add_model_specific_part(const Json_Site &other) = 0;

    public:
        Json_Site() {
            for (const auto &element : gram::json::spec::site_fields.items()) {
                json_site.emplace(element.key(), JSON::array());
            }
        }

        Json_Site(Json_Site const &other) : json_site(other.json_site) {}

        // Functions implementing site combining
        void build_allele_combi_map(JSON const &json_site, allele_combi_map &m);

        void append_entries_from(JSON const &json_site);

        allele_vec get_all_alleles(allele_combi_map &m);

        JSON rescale_entries(allele_combi_map const &m) const;

        void combine_with(const Json_Site &other);

        JSON const &get_site() const { return json_site; }

        JSON get_site_copy() const { return json_site; }

        void set_site(JSON const &json_site) { this->json_site = json_site; }
    };

    class LevelGenotyped_Json_Site : public Json_Site {
    public:
        LevelGenotyped_Json_Site() : Json_Site() {
            json_site.emplace("GT_CONF", JSON::array());
        }

        void add_model_specific_part(const Json_Site &other) override {};
    };
}

#endif //SITE_JSON_SPEC
