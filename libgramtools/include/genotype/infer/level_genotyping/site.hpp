#ifndef LVLGT_SITE
#define LVLGT_SITE

#include <variant>
#include "genotype/infer/interfaces.hpp"
#include "genotype/infer/output_specs/json_site_spec.hpp"

namespace gram::genotype::infer {

    class LevelGenotypedSite : public GenotypedSite {
        double gt_conf = 0.; /**< Difference in log likelihood between most likely and next most likely genotype **/
        double gt_conf_percentile = 0.; /**< Percent of gt_confs that are below this value */
    public:
        ~LevelGenotypedSite() override = default;

        void set_gt_conf(double const& val) {this->gt_conf = val;}
        void set_gt_conf_percentile(double const& val) {this->gt_conf_percentile = val;}

        site_entries get_model_specific_entries() override;
        void null_model_specific_entries() override;
    };
}

#endif //LVLGT_SITE
