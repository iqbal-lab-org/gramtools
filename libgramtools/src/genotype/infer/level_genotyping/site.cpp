#include "genotype/infer/level_genotyping/site.hpp"

using namespace gram::genotype::infer;

JSON LevelGenotypedSite::get_JSON() {
    if (! site_json.at("GT").empty()) return site_json;

    for (int i{0}; i < alleles.size(); ++i) site_json.at("ALS").push_back(alleles.at(i).sequence);

    if (is_null()) site_json.at("GT").push_back(JSON::array({nullptr}));
    else site_json.at("GT").push_back(JSON::array({std::get<GtypedIndices>(genotype)}));

    site_json.at("HAPG").push_back(JSON::array({haplogroups}));

    site_json.at("COVS").push_back(JSON::array({allele_covs}));

    site_json.at("DP").push_back(total_coverage);

    site_json["GT_CONF"] = JSON::array({gt_conf});
    return site_json;
}
