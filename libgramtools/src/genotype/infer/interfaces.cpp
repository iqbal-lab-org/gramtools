#include "genotype/infer/interfaces.hpp"
#include "genotype/infer/output_specs/json_prg_spec.hpp"
#include "genotype/infer/output_specs/json_site_spec.hpp"

namespace gram::genotype::infer {


    void GenotypedSite::populate_site(gtype_information const& gtype_info){
        this->gtype_info.alleles = gtype_info.alleles;
        this->gtype_info.genotype = gtype_info.genotype;
        this->gtype_info.allele_covs = gtype_info.allele_covs;
        this->gtype_info.total_coverage = gtype_info.total_coverage;
        this->gtype_info.haplogroups = gtype_info.haplogroups;
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
    assert(gtype_info.alleles.size() > 0);
    assert(num_haplogroups > 0);
    AlleleIds result;

    AlleleIdSet genotyped_haplogroups;
    for (auto const& gt : std::get<GtypedIndices>(gtype_info.genotype)){
       genotyped_haplogroups.insert(gtype_info.alleles.at(gt).haplogroup);
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

void Genotyper::run_invalidation_process(gt_site_ptr const& genotyped_site, Marker const& site_ID) {
    // Invalidation process attempted only if this site contains 1+ site
    if (!genotyped_site->is_null() && child_m.find(site_ID) != child_m.end()) {
        auto candidate_haplogroups =
                genotyped_site->get_nonGenotyped_haplogroups();
        auto haplogroups_with_sites = get_haplogroups_with_sites(site_ID, candidate_haplogroups);
        invalidate_if_needed(site_ID, haplogroups_with_sites);
    }
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

    for (int i{0}; i < gtype_info.alleles.size(); ++i) json_site_copy.at("ALS")
        .push_back(gtype_info.alleles.at(i).sequence);

    if (is_null()) json_site_copy.at("GT").push_back(JSON::array({nullptr}));
    else json_site_copy.at("GT").push_back(JSON(std::get<GtypedIndices>(gtype_info.genotype)));

    json_site_copy.at("HAPG").push_back(JSON(gtype_info.haplogroups));

    json_site_copy.at("COV").push_back(JSON(gtype_info.allele_covs));
    json_site_copy.at("DP").push_back(gtype_info.total_coverage);

    this->add_model_specific_JSON(json_site_copy);
    json_site->set_site(json_site_copy);

    return json_site;
}
}
