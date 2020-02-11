#include "genotype/infer/level_genotyping/site.hpp"

using namespace gram::genotype::infer;

void LevelGenotypedSite::add_JSON() {
    site_json.at("COVS").push_back(JSON(allele_covs));

    site_json.at("DP").push_back(total_coverage);

    site_json["GT_CONF"] = JSON::array({gt_conf});
}
