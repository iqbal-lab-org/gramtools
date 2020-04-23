#include "genotype/infer/level_genotyping/site.hpp"
#include "genotype/infer/output_specs/coverage_common.hpp"

using namespace gram::genotype::infer;

site_entries LevelGenotypedSite::get_model_specific_entries(){
    site_entry<double> gt_conf_entry {
           "FORMAT", "GT_CONF", {gt_conf}, true
            };
    site_entries result{{gt_conf_entry}};
    return result;
}

void LevelGenotypedSite::null_model_specific_entries() {
    this->gt_conf = 0.;
}
