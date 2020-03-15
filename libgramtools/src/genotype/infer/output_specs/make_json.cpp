#include "genotype/infer/output_specs/make_json.hpp"
#include "prg/coverage_graph.hpp"

json_prg_ptr make_json_prg(gtyper_ptr const& gtyper){
    auto result = std::make_shared<Json_Prg>();
    populate_json_prg(*result, gtyper);
    auto genotyped_records = gtyper->get_genotyped_records();
    for (auto const& site : genotyped_records)
        result->add_site(make_json_site(site));
    return result;
}

void populate_json_prg(Json_Prg& json_prg, gtyper_ptr const& gtyper){
    JSON json_prg_copy = json_prg.get_prg_copy();
    auto cov_graph = gtyper->get_cov_g();
    auto child_m = gtyper->get_child_m();
    auto genotyped_records = gtyper->get_genotyped_records();
    if (! cov_graph->is_nested) json_prg_copy.at("Lvl1_Sites").push_back("all");
    else {
        for (int i{0}; i < genotyped_records.size(); ++i)
            if (cov_graph->par_map.find(index_to_siteID(i)) ==
                cov_graph->par_map.end())
                json_prg_copy.at("Lvl1_Sites").push_back(i);

        for (const auto &child_entry : child_m) {
            auto site_index = std::to_string(siteID_to_index(child_entry.first));
            json_prg_copy.at("Child_Map").emplace(site_index, JSON::object());
            for (const auto &hapg_entry : child_entry.second) {
                auto copy = hapg_entry.second;
                for (auto &el : copy) el = siteID_to_index(el);
                json_prg_copy.at("Child_Map").at(site_index)[std::to_string(hapg_entry.first)] =
                        JSON(copy);
            }
        }
    }
    json_prg.set_prg(json_prg_copy);

    for (auto const& header : gtyper->get_model_specific_headers())
        json_prg.add_header(header);
}

void add_model_specific_entries(JSON& json_site, site_entries const& entries){
   for (auto const& entry: entries.doubles){
       json_site[entry.ID] = JSON::array();
      if (entry.single_val) json_site.at(entry.ID).push_back(entry.vals.at(0));
      else json_site.at(entry.ID).push_back(JSON(entry.vals));
   }
}


json_site_ptr make_json_site(gt_site_ptr const& gt_site){
    json_site_ptr result = std::make_shared<Json_Site>();
    auto json_site = result->get_site_copy();
    auto gtype_info = gt_site->get_all_gtype_info();
    for (int i{0}; i < gtype_info.alleles.size(); ++i) json_site.at("ALS")
                .push_back(gtype_info.alleles.at(i).sequence);

    if (gt_site->is_null()) json_site.at("GT").push_back(JSON::array({nullptr}));
    else json_site.at("GT").push_back(JSON(gtype_info.genotype));

    json_site.at("HAPG").push_back(JSON(gtype_info.haplogroups));

    json_site.at("COV").push_back(JSON(gtype_info.allele_covs));
    json_site.at("DP").push_back(gtype_info.total_coverage);

    add_model_specific_entries(json_site, gt_site->get_model_specific_entries());
    result->set_site(json_site);

    return result;
}

