#include "genotype/infer/output_specs/json_prg_spec.hpp"
#include "genotype/infer/output_specs/json_site_spec.hpp"

using namespace gram::json;

Json_Prg::Json_Prg(JSON input_json) : json_prg(std::move(input_json)) {
  for (JSON const& site_json : json_prg.at("Sites")) {
    sites.emplace_back(std::make_shared<Json_Site>(site_json));
  }
}

void Json_Prg::set_prg(JSON const& input_json) {
  json_prg = input_json;
  for (JSON const& site_json : json_prg.at("Sites")) {
    sites.emplace_back(std::make_shared<Json_Site>(site_json));
  }
}

void Json_Prg::set_sample_info(std::string const& name,
                               std::string const& desc) {
  if (json_prg.at("Samples").size() > 1)
    throw JSONConsistencyException("This JSON already contains > 1 samples");

  json_prg.at("Samples") = JSON::array();
  json_prg.at("Samples").push_back(JSON{{"Name", name}, {"Desc", desc}});
}

void Json_Prg::add_site(json_site_ptr const& json_site) {
  sites.push_back(json_site);
  json_prg.at("Sites").push_back(json_site->get_site());
}

void Json_Prg::add_header(vcf_meta_info_line header) {
  if (header.meta_type == "Model") {
    json_prg.at("Model") = header.flat_value;
  } else if (header.meta_type == "FORMAT") {
    json_prg.at("Site_Fields")[header.ID] = {{"Desc", header.desc}};
  }
}

void Json_Prg::add_samples(Json_Prg& other, const bool force) {
  auto& other_prg = other.get_prg();
  if (other_prg.at("Sites").at(0).at("GT").size() !=
      other_prg.at("Samples").size())
    throw JSONConsistencyException(
        "Merged in JSON does not have number of GT arrays"
        " consistent with its number of Samples");

  std::map<std::string, std::size_t> duplicates;
  for (auto const& e : json_prg.at("Samples"))
    duplicates.insert({e.at("Name"), 1});

  for (auto& sample_entry : other_prg.at("Samples")) {
    std::string const name = sample_entry.at("Name");
    std::string used_name = name;
    if (duplicates.find(name) != duplicates.end()) {
      if (!force)
        throw JSONConsistencyException(
            std::string{"Duplicate sample name found: " + name});
      else {
        auto& time_seen = duplicates.at(name);
        used_name = name + std::string("_") + std::to_string(time_seen);
        time_seen++;
      }
    } else
      duplicates.insert({name, 1});

    sample_entry.at("Name") = used_name;
    json_prg.at("Samples").push_back(sample_entry);
  }
}

void Json_Prg::combine_with(Json_Prg& other, bool force) {
  auto other_prg = other.get_prg();
  if (json_prg.at("Model") != other_prg.at("Model"))
    throw JSONCombineException("JSONs have different models");

  auto bad_prg = json_prg.at("Lvl1_Sites") != other_prg.at("Lvl1_Sites");
  bad_prg |= json_prg.at("Child_Map") != other_prg.at("Child_Map");

  if (bad_prg)
    throw JSONCombineException(
        "Incompatible PRGs (Check Child_Map and Lvl1_Sites)");

  if (json_prg.at("Site_Fields") != other_prg.at("Site_Fields"))
    throw JSONCombineException("Incompatible Site Fields");

  if (sites.size() != other.sites.size())
    throw JSONCombineException("JSONs do not have the same number of sites");

  add_samples(other, force);

  auto const gtyping_model = json_prg.at("Model");
  for (std::size_t j{0}; j < sites.size(); j++) {
    sites.at(j)->combine_with(*other.sites.at(j), gtyping_model);
    json_prg.at("Sites").at(j) = sites.at(j)->get_site();
  }
}
