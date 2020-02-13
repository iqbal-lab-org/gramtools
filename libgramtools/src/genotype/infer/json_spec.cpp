#include "genotype/infer/json_spec.hpp"

using namespace gram::json;

void Json_Prg::add_site(json_site_ptr json_site){
   sites.push_back(json_site);
   json_prg.at("Sites").push_back(json_site->get_site_copy());
}

void Json_Prg::combine_with(const Json_Prg& other){
    auto other_prg = other.get_prg();
    if (json_prg.at("Model") != other_prg.at("Model"))
        throw JSONCombineException("JSONs have different models");

    auto bad_prg = json_prg.at("Lvl1_Sites") != other_prg.at("Lvl1_Sites");
    bad_prg |= json_prg.at("Child_Map") != other_prg.at("Child_Map");

    if (bad_prg) throw JSONCombineException("Incompatible PRGs (Check Child_Map and Lvl1_Sites)");

    if (json_prg.at("Site_Fields") != other_prg.at("Site_Fields"))
        throw JSONCombineException("Incompatible Site Fields");

    if (sites.size() != other.sites.size())
        throw JSONCombineException("JSONs do not have the same number of sites");
}

void Json_Site::combine_with(const Json_Site &other){
    auto other_site = other.get_site_copy();
    if (json_site.at("ALS").at(0) != other_site.at("ALS").at(0)){
        std::string msg("Sites do not have same 'reference' allele: ");
        msg = msg + std::string(json_site.at("ALS").at(0)) + std::string(" vs ") +
                std::string(other_site.at("ALS").at(0));
        throw JSONCombineException(msg);
    }

    allele_combi_map m;
    build_allele_combi_map(json_site, m);
    build_allele_combi_map(other_site, m);
}


void Json_Site::build_allele_combi_map(JSON const& json_site, allele_combi_map& m){
    using namespace gram::genotype::infer;

    auto insertion_index{m.size()};
    auto const num_samples{json_site.at("GT").size()};

    std::size_t sample_num{0};
    while (sample_num < num_samples){
        GtypedIndices const& gts = json_site.at("GT").at(sample_num);
        AlleleIds const& hapgs = json_site.at("HAPG").at(sample_num);
        if (gts.size() != hapgs.size())
            throw JSONConsistencyException("Different number of GT and HAPG entries");

        for (std::size_t j{0}; j < gts.size(); j++) {
           auto const hapg = hapgs.at(j);
           auto const allele = json_site.at("ALS").at(gts.at(j));
           if (m.find(allele) == m.end()){
               m.insert({allele, site_rescaler{insertion_index++, hapg}});
           }
           else if (m.at(allele).hapg != hapg){
                   std::string msg("Allele has two HAPG values: ");
                   msg = msg + std::to_string(hapg) +
                           std::string(" vs ") +
                           std::to_string(m.at(allele).hapg);
                   throw JSONConsistencyException(msg);
               }
        }
        sample_num++;
    }
}
