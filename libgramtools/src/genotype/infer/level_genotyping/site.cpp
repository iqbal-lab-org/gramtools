#include "genotype/infer/level_genotyping/site.hpp"

using namespace gram::genotype::infer;

void LevelGenotypedSite::add_model_specific_JSON(JSON& input_json) {
    input_json["GT_CONF"] = JSON::array({gt_conf});
}
