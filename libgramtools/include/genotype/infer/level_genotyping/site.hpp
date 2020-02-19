#ifndef LVLGT_SITE
#define LVLGT_SITE

#include <variant>
#include "common/utils.hpp"
#include "genotype/infer/interfaces.hpp"
#include "genotype/infer/json_spec/site_spec.hpp"

namespace gram::genotype::infer {

    class LevelGenotypedSite : public GenotypedSite {
        double gt_conf = 0.; /**< Difference in log likelihood between most likely and next most likely genotype **/
    public:
        LevelGenotypedSite() {
            auto json_site_ptr = std::make_shared<LevelGenotyped_Json_Site>();
            json_site = json_site_ptr;
        }

        ~LevelGenotypedSite() override = default;


        void set_gt_conf(double const& gt_conf) {this->gt_conf = gt_conf;}
        void add_model_specific_JSON(JSON& input_json) override;
    };
}

#endif //LVLGT_SITE
