#include "genotype/infer/level_genotyping/site.hpp"
#include "genotype/infer/output_specs/json_common.hpp"

using namespace gram::genotype::infer;

entry_vec LevelGenotypedSite::get_model_specific_entries(){
    entry_vec result;
    auto gt_conf_entry = std::make_shared<base_site_entry>(site_entry<double>{
           "FORMAT", "GT_CONF", {gt_conf}, true
            });
    result.push_back(gt_conf_entry);
    return result;
//    input_json["GT_CONF"] = JSON::array({gt_conf});
}

void LevelGenotypedSite::null_model_specific_entries() {
    this->gt_conf = 0.;
}
