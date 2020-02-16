#include "genotype/infer/interfaces.hpp"
#include "genotype/infer/json_spec/prg_spec.hpp"
#include "genotype/infer/json_spec/site_spec.hpp"

namespace gram::genotype::infer {

void Genotyper::populate_json_prg() {
    JSON json_prg_copy = json_prg->get_prg_copy();
    if (! cov_graph->is_nested) json_prg_copy.at("Lvl1_Sites").push_back("all");
    else {
        for (int i{0}; i < genotyped_records.size(); ++i)
            if (cov_graph->par_map.find(index_to_siteID(i)) ==
                cov_graph->par_map.end())
                json_prg_copy.at("Lvl1_Sites").push_back(i);

        for (const auto &child_entry : child_m) {
            auto site_index = std::to_string(siteID_to_index(child_entry.first));
            json_prg_copy.at("Child_map").emplace(site_index, JSON::object());
            for (const auto &hapg_entry : child_entry.second) {
                auto copy = hapg_entry.second;
                for (auto &el : copy) el = siteID_to_index(el);
                json_prg_copy.at("Child_map").at(site_index)[std::to_string(hapg_entry.first)] =
                        JSON(copy);
            }
        }
    }
    json_prg->set_prg(json_prg_copy);
}

void Genotyper::add_json_sites(){
    for (auto const& site : genotyped_records)
        json_prg->add_site(site->get_JSON());
}

    JSON Genotyper::get_JSON() {
    if (json_prg->get_prg().at("Sites").empty()){
        populate_json_prg();
        add_json_sites();
    }
    return json_prg->get_prg();
}

allele_vector const GenotypedSite::get_unique_genotyped_alleles(allele_vector const &all_alleles,
                                                                GenotypeOrNull const &genotype) const {

    std::set<GtypedIndex> distinct_genotypes;
    if (auto valid_gtype = std::get_if<GtypedIndices>(&genotype)) {
        // NOTE/CRUCIAL: this sorts the genotypes (eg 1,0 goes to 0, 1), which is REQUIRED for REF allele production
        distinct_genotypes = std::set<GtypedIndex>(valid_gtype->begin(), valid_gtype->end());
    } else distinct_genotypes = std::set<GtypedIndex>{0}; // If null genotype, take the reference only

    allele_vector result(distinct_genotypes.size());

    std::size_t insertion_index{0};
    for (auto const &allele_id : distinct_genotypes) {
        result.at(insertion_index) = all_alleles.at(allele_id);
        ++insertion_index;
    }
    return result;
}

AlleleIds const GenotypedSite::get_nonGenotyped_haplogroups() const{
    assert(! is_null());
    assert(alleles.size() > 0);
    assert(num_haplogroups > 0);
    AlleleIds result;

    AlleleIdSet genotyped_haplogroups;
    for (auto const& gt : std::get<GtypedIndices>(genotype)){
       genotyped_haplogroups.insert(alleles.at(gt).haplogroup);
    }

    for (AlleleId i{0}; i < num_haplogroups; i++){
        if (genotyped_haplogroups.find(i) == genotyped_haplogroups.end())
            result.push_back(i);
    }
    return result;
}

AlleleIds GenotypedSite::get_genotyped_haplogroups(allele_vector const& input_alleles,
                                                   GtypedIndices const& input_gts) const{
    AlleleIds result;
    for (auto const& gt : input_gts){
       result.push_back(input_alleles.at(gt).haplogroup);
    }
    return result;
}

AlleleIds Genotyper::get_haplogroups_with_sites(Marker const& site_ID, AlleleIds candidate_haplogroups) const{
    AlleleIds result{};
    if (child_m.find(site_ID) == child_m.end()) return result;
    auto child_entry = child_m.at(site_ID);
    for (auto const& candidate : candidate_haplogroups){
        if (child_entry.find(candidate) != child_entry.end()) result.push_back(candidate);
    }
    return result;
}

void Genotyper::invalidate_if_needed(Marker const& parent_site_ID, AlleleIds haplogroups){
    if (haplogroups.size() == 0) return;

    std::vector<VariantLocus> to_process;
    for (auto const& haplogroup : haplogroups) to_process.push_back(VariantLocus{parent_site_ID, haplogroup});

    VariantLocus cur_locus;
    while (to_process.size() > 0){
        cur_locus = to_process.back();
        to_process.pop_back();

        // We know the site/haplogroup combination bears 1+ site so use .at() .
        auto sites_on_haplogroup = child_m.at(cur_locus.first).at(cur_locus.second);
        for (auto const& child_marker : sites_on_haplogroup){
            auto referent_site = genotyped_records.at(siteID_to_index(child_marker));
            if (referent_site->is_null()) continue;
            referent_site->make_null();

            auto all_haplos = referent_site->get_all_haplogroups();
            auto haplos_with_sites = get_haplogroups_with_sites(child_marker, all_haplos);
            for (auto const& h : haplos_with_sites) to_process.push_back(VariantLocus{child_marker,h});
        }
    }
}

json_site_ptr GenotypedSite::get_JSON(){
    auto json_site_copy = json_site->get_site_copy();
    if (! json_site_copy.at("GT").empty()) return json_site;

    for (int i{0}; i < alleles.size(); ++i) json_site_copy.at("ALS").push_back(alleles.at(i).sequence);

    if (is_null()) json_site_copy.at("GT").push_back(JSON::array({nullptr}));
    else json_site_copy.at("GT").push_back(JSON(std::get<GtypedIndices>(genotype)));

    json_site_copy.at("HAPG").push_back(JSON(haplogroups));

    json_site_copy.at("COVS").push_back(JSON(allele_covs));
    json_site_copy.at("DP").push_back(total_coverage);

    this->add_model_specific_JSON(json_site_copy);
    json_site->set_site(json_site_copy);

    return json_site;
}
}
