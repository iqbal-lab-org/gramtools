#include "genotype/infer/json_spec.hpp"

using namespace gram::json;

void Json_Prg::add_site(json_site_ptr json_site){
   sites.push_back(json_site);
   json_prg.at("Sites").push_back(json_site->get_site_copy());
}

