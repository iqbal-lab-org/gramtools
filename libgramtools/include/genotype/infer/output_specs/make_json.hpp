#ifndef GTYPE_MAKE_JSON_HPP
#define GTYPE_MAKE_JSON_HPP

#include "json_prg_spec.hpp"
#include "json_site_spec.hpp"
#include "genotype/infer/interfaces.hpp"

using namespace gram::json;
using namespace gram::json::spec;
using namespace gram::genotype::infer;

json_site_ptr make_json_site(gt_site_ptr const& gt_site);

/**
 * Populates the PRG-related entries (Lvl1_sites, child map) of a Json_Prg class.
 */
void populate_json_prg(Json_Prg const& json_prg, gtyper_ptr const& gtyper);

Json_Prg make_json_prg(gtyper_ptr const& gtyper);


#endif //GTYPE_MAKE_JSON_HPP
