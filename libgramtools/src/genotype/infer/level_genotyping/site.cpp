#include "genotype/infer/level_genotyping/site.hpp"

using namespace gram::genotype::infer;

JSON LevelGenotypedSite::get_JSON() {
    if (! site_json["GT"].empty()) return site_json;
    JSON result;

    for (int i{0}; i < alleles.size(); ++i) result["alleles"].push_back(alleles.at(i).sequence);

    if (is_null()) result["GT"].push_back(JSON::array({nullptr}));
    else result["GT"].push_back(JSON::array({std::get<GtypedIndices>(genotype)}));

    result["GT_CONF"] = JSON::array({gt_conf});

    result["HAPG"] = JSON::array();
    result["HAPG"].push_back(JSON::array({haplogroups}));

    result["COVS"] = JSON::array();
    result["COVS"].push_back(JSON::array({allele_covs}));

    result["DP"] = JSON::array({total_coverage});

    site_json = result;
    return result;
}
