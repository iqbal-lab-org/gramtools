#ifndef SITE_JSON_SPEC
#define SITE_JSON_SPEC

#include "common/data_types.hpp"
#include "fields.hpp"

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
    auto site_fields = spec::site_fields();
    for (const auto &element : site_fields.items()) {
      json_site.emplace(element.key(), JSON::array());
    }
    json_site.emplace("SEG", "");
  }

  Json_Site(JSON const input_json) : json_site(input_json) {}

  // ____Functions implementing site combining____
  void build_allele_combi_map(JSON const &json_site, allele_combi_map &m);
  allele_vec get_all_alleles(allele_combi_map &m);
  void rescale_entries(allele_combi_map const &m);
  void append_trivial_entries_from(JSON const &input_site);
  void add_model_specific_entries_from(JSON const &input_site,
                                       std::string const gtyping_model);
  void combine_with(Json_Site &other, std::string const &gtyping_model = "");

  JSON const &get_site() const { return json_site; }
  JSON &get_site() { return json_site; }

  void set_site(JSON const input_site) { json_site = input_site; }
  void set_pos(std::size_t pos) { json_site.at("POS") = pos; }
  void set_segment(std::string ID) { this->json_site.at("SEG") = ID; }
};

}  // namespace gram::json

#endif  // SITE_JSON_SPEC
