/** @file
 * vcf and JSON/jvcf output specs
 */

#ifndef GTYPE_FIELDS_HPP
#define GTYPE_FIELDS_HPP

#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

using strings = std::vector<std::string>;

namespace gram::genotype::output_spec {

struct vcf_meta_info_line {
  std::string meta_type, ID = "", desc = "", flat_value = "", num = "",
                         type = "";
  std::size_t length{0};

  /**
   * Constructor for generic header, simple key-val. Eg: ##source=my_source
   */
  vcf_meta_info_line(std::string m_t, std::string flat_val)
      : meta_type(m_t), flat_value(flat_val) {}

  /**
   * Constructor for structured header of type: FORMAT, INFO
   */
  vcf_meta_info_line(std::string m_t, std::string id, std::string desc,
                     std::string num, std::string type)
      : meta_type(m_t), ID(id), desc(desc), num(num), type(type) {}

  /**
   * Constructor for structured header of type: FILTER, ALT
   */
  vcf_meta_info_line(std::string m_t, std::string id, std::string desc)
      : meta_type(m_t), ID(id), desc(desc) {}

  /**
   * Constructor for structured header of type: contig
   */
  vcf_meta_info_line(std::string m_t, std::string id, std::size_t len)
      : meta_type(m_t), ID(id), length(len) {}

  std::string const to_string() const {
    std::stringstream out;
    out << "##" << meta_type << "=";
    if (flat_value != "") {
      out << flat_value;
      return out.str();
    }
    out << "<ID=" << ID;
    if (num != "") out << ",Number=" << num;
    if (type != "") out << ",Type=" << type;
    if (desc != "") out << ",Description=\"" << desc << "\"";
    if (length != 0) out << ",length=" << std::to_string(length);
    out << ",Source=\"gramtools\"";
    out << ">";
    return out.str();
  }
};

using header_vec = std::vector<vcf_meta_info_line>;

template <typename T>
struct site_entry {
  site_entry() = default;
  site_entry(vcf_meta_info_line const& spec)
      : meta_type(spec.meta_type), ID(spec.ID) {}
  std::string meta_type, ID;
  std::vector<T> vals;
  bool single_val = true;
};

/*
 * I cannot figure out a nice way to store a set of arbitrarily typed site_entry
 * objects, so defaulted to explicitly using those that have come up
 */
struct site_entries {
  std::vector<site_entry<double>> doubles = {};
};

static header_vec const common_headers{
    vcf_meta_info_line{"FORMAT", "GT", "Genotype", "1", "String"},
    vcf_meta_info_line{"FORMAT", "DP", "Total read depth on variant site", "1",
                       "Integer"},
    vcf_meta_info_line{"FORMAT", "COV", "Read coverage on each allele", "R",
                       "Float"},
    vcf_meta_info_line{"FORMAT", "FT", "Filters failed in a sample", "1",
                       "String"},
    vcf_meta_info_line{"FILTER", "AMBIG",
                       "Ambiguous site. Different variant paths can produce "
                       "the same sequence."}};
}  // namespace gram::genotype::output_spec

using JSON = nlohmann::json;
using namespace gram::genotype::output_spec;

namespace gram::json {
class Json_Prg;
class Json_Site;

static const strings trivially_merged_entries{"GT", "HAPG", "COV", "DP", "FT"};
static const strings singleton_entries{"POS", "SEG"};

using json_prg_ptr = std::shared_ptr<Json_Prg>;
using json_site_ptr = std::shared_ptr<Json_Site>;
using json_site_vec = std::vector<json_site_ptr>;

class JSONCombineException : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

class JSONConsistencyException : public std::runtime_error {
  using std::runtime_error::runtime_error;
};
}  // namespace gram::json

namespace gram::json::spec {

static JSON site_fields() {
  JSON result = {
      {"POS", {{"Desc", "Position on reference or pseudo-reference"}}},
      {"SEG", {{"Desc", "Segment ID"}}},
      {"ALS", {{"Desc", "Alleles at this site"}}},
      {
          "HAPG",
          {{"Desc", "Sample haplogroups of genotyped alleles"}},
      },
  };

  // Populate with headers common with vcf output
  for (auto const& header : common_headers) {
    if (header.meta_type == "FORMAT")
      result[header.ID] = {{"Desc", header.desc}};
  }
  return result;
}

static JSON filters() {
  JSON result;
  // Populate with headers common with vcf output
  for (auto const& header : common_headers) {
    if (header.meta_type == "FILTER")
      result[header.ID] = {{"Desc", header.desc}};
  }
  return result;
}

const JSON json_prg{
    {"Model", "UNKNOWN"},         {"Site_Fields", site_fields()},
    {"Filters", filters()},       {"Samples", JSON::array()},
    {"Sites", JSON::array()},     {"Lvl1_Sites", JSON::array()},
    {"Child_Map", JSON::object()}};
}  // namespace gram::json::spec

#endif  // GTYPE_FIELDS_HPP
