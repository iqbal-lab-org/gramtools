#include "genotype/infer/output_specs/make_json.hpp"
#include "prg/coverage_graph.hpp"

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
}

Json_Prg make_json_prg(gtyper_ptr const& gtyper){
    Json_Prg json_prg;
    populate_json_prg(json_prg, gtyper);
    auto genotyped_records = gtyper->get_genotyped_records();
    for (auto const& site : genotyped_records)
        json_prg.add_site(site->get_JSON());
    return json_prg;
}

