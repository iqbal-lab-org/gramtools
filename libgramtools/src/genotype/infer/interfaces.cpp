#include "genotype/infer/interfaces.hpp"

namespace gram::genotype::infer {

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

JSON GenotypedSite::get_JSON(){
    if (! site_json.at("GT").empty()) return site_json;

    for (int i{0}; i < alleles.size(); ++i) site_json.at("ALS").push_back(alleles.at(i).sequence);

    if (is_null()) site_json.at("GT").push_back(JSON::array({nullptr}));
    else site_json.at("GT").push_back(JSON(std::get<GtypedIndices>(genotype)));

    site_json.at("HAPG").push_back(JSON(haplogroups));

    this->add_JSON();

    return site_json;
}
}
