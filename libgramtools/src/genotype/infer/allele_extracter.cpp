#include "genotype/infer/allele_extracter.hpp"
#include "genotype/infer/interfaces.hpp"
#include "prg/coverage_graph.hpp"

using namespace gram::genotype::infer;

AlleleExtracter::AlleleExtracter(covG_ptr site_start, covG_ptr site_end,
                                 gt_sites& sites)
    : genotyped_sites(&sites) {
  assert(site_start->is_bubble_start());
  AlleleId haplogroup_ID{FIRST_ALLELE};

  for (auto& haplogroup_start_node : site_start->get_edges()) {
    allele_vector extracted_alleles =
        extract_alleles(haplogroup_ID, haplogroup_start_node, site_end);
    alleles.insert(alleles.end(), extracted_alleles.begin(),
                   extracted_alleles.end());
    haplogroup_ID++;
  }
}

allele_vector AlleleExtracter::allele_combine(allele_vector const& existing,
                                              std::size_t site_index) {
  // Sanity check: site_index refers to actual site
  assert(0 <= site_index && site_index < genotyped_sites->size());
  gt_site_ptr referent_site = genotyped_sites->at(site_index);

  allele_vector relevant_alleles =
      referent_site->get_unique_genotyped_alleles();
  if (referent_site->extra_alleles()) {
    allele_vector const& to_add = referent_site->extra_alleles().value();
    relevant_alleles.insert(relevant_alleles.end(), to_add.begin(),
                            to_add.end());
  }
  if (relevant_alleles.empty())
    relevant_alleles.push_back(referent_site->get_alleles().at(0));
  allele_vector combinations(existing.size() * relevant_alleles.size());

  std::size_t insertion_index{0};
  for (auto const& allele : existing) {
    for (auto const& added_allele : relevant_alleles) {
      combinations.at(insertion_index) = allele + added_allele;
      ++insertion_index;
    }
  }
  return combinations;
}

void AlleleExtracter::allele_paste(allele_vector& existing,
                                   covG_ptr sequence_node) {
  Allele to_paste_allele{sequence_node->get_sequence(),
                         sequence_node->get_coverage()};
  for (auto& allele : existing) allele = allele + to_paste_allele;
}

void AlleleExtracter::place_ref_as_first_allele(allele_vector& alleles,
                                                Allele const& ref_allele) {
  this->_ref_allele_got_made_naturally = false;
  auto found_ref = std::find(alleles.begin(), alleles.end(), ref_allele);
  if (found_ref == alleles.end()) {
    alleles = prepend(alleles, ref_allele);
  } else if (found_ref != alleles.begin())
    std::swap(*found_ref, alleles.at(0));
  else
    this->_ref_allele_got_made_naturally = true;
}

allele_vector AlleleExtracter::extract_alleles(AlleleId const haplogroup,
                                               covG_ptr haplogroup_start,
                                               covG_ptr site_end) {
  allele_vector haplogroup_alleles{
      {"", {}, haplogroup}};  // Make one empty allele as starting point, allows
                              // for direct deletion
  covG_ptr cur_Node{haplogroup_start};

  Allele ref_allele{"", {}, 0};

  while (cur_Node != site_end) {
    if (cur_Node->is_bubble_start()) {
      auto site_index = siteID_to_index(cur_Node->get_site_ID());
      haplogroup_alleles = allele_combine(haplogroup_alleles, site_index);

      auto referent_site = genotyped_sites->at(site_index);
      cur_Node =
          referent_site->get_site_end_node();  // Move past site, to bubble end

      if (haplogroup == 0) {  // REF extraction
        auto ref_containing_alleles = referent_site->get_alleles();
        assert(ref_containing_alleles.size() > 0);
        ref_allele = ref_allele + ref_containing_alleles.at(0);
      }
    } else {
      allele_paste(haplogroup_alleles, cur_Node);
      if (haplogroup == 0)
        ref_allele = ref_allele +
                     Allele{cur_Node->get_sequence(), cur_Node->get_coverage()};
    }

    // The only nodes with >1 neighbour are bubble starts and we
    // skipped past those.
    assert(cur_Node->get_num_edges() == 1);

    cur_Node = *(cur_Node->get_edges().begin());  // Advance to the next node
  }

  if (haplogroup == 0)
    place_ref_as_first_allele(haplogroup_alleles, ref_allele);

  return haplogroup_alleles;
}
