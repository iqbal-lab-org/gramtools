#include "genotype/infer/json_spec/site_spec.hpp"

using namespace gram::json;
using namespace gram::genotype::infer;

void add_or_check_allele(std::string const& allele, gram::AlleleId const& hapg, allele_combi_map& m,
                         std::size_t& insertion_index){
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

void Json_Site::build_allele_combi_map(JSON const& json_site, allele_combi_map& m){

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
            add_or_check_allele(allele, hapg, m, insertion_index);
        }
        sample_num++;
    }
}

allele_vec Json_Site::get_all_alleles(allele_combi_map& m){
    allele_vec result(m.size());
    for (auto const& entry : m){
        result.at(entry.second.index) = entry.first;
    }
    return result;
}

JSON Json_Site::rescale_entries(allele_combi_map const &m) const{
    auto result = json_site;
    auto const num_samples{result.at("GT").size()};

    std::size_t sample_num{0};
    while (sample_num < num_samples){
        GtypedIndices gts = result.at("GT").at(sample_num);

        allele_coverages const covs = json_site.at("COVS").at(sample_num);
        allele_coverages new_covs(m.size(), 0);

        if (json_site.at("ALS").size() != covs.size())
            throw JSONConsistencyException("Different number of ALS and COVS entries");

        for (auto& gt : gts) {
            auto const allele = json_site.at("ALS").at(gt);
            auto const rescaled_idx = m.at(allele).index;
            gt = rescaled_idx;
        }
        for (int j{0}; j < covs.size(); j++){
            auto const allele = json_site.at("ALS").at(j);
            // The allele is not called in any sample, so we ignore it
            if (m.find(allele) == m.end()) continue;
            auto const rescaled_idx = m.at(allele).index;
            new_covs.at(rescaled_idx) = covs.at(j);
        }
        result.at("GT").at(sample_num) = gts;
        result.at("COVS").at(sample_num) = new_covs;
        sample_num++;
    }
    return result;
}

void Json_Site::append_entries_from(JSON const& json_site){
    std::vector<std::string> entries{"GT", "HAPG", "COVS", "DP"};
    for (auto const& entry : entries){
        for (auto const& element : json_site.at(entry))
            this->json_site.at(entry).push_back(element);
    }
}

void Json_Site::combine_with(const Json_Site &other){
    auto other_json = other.get_site_copy();
    std::string this_ref = json_site.at("ALS").at(0);
    if (this_ref != other_json.at("ALS").at(0)){
        std::string msg("Sites do not have same 'reference' allele: ");
        msg = msg + std::string(json_site.at("ALS").at(0)) + std::string(" vs ") +
              std::string(other_json.at("ALS").at(0));
        throw JSONCombineException(msg);
    }

    allele_combi_map m{{this_ref, site_rescaler{0, 0}}}; // Always place the REF
    build_allele_combi_map(json_site, m);
    build_allele_combi_map(other_json, m);

    auto this_rescaled = rescale_entries(m);
    this_rescaled.at("ALS") = get_all_alleles(m);
    set_site(this_rescaled);
    JSON other_rescaled = other.rescale_entries(m);

    append_entries_from(other_rescaled);
}

