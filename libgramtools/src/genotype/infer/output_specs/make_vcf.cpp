#include "genotype/infer/output_specs/make_vcf.hpp"
#include <htslib/synced_bcf_reader.h>
#include "genotype/infer/output_specs/fields.hpp"
#include "genotype/infer/output_specs/segment_tracker.hpp"
#include "prg/coverage_graph.hpp"

void write_vcf(gram::GenotypeParams const& params, gtyper_ptr const& gtyper,
               SegmentTracker& tracker) {
  auto fout = bcf_open(params.genotyped_vcf_fpath.c_str(), "wz");  // Writer

  // Set up and write header
  bcf_hdr_t* header = bcf_hdr_init("w");
  populate_vcf_hdr(header, gtyper, params, tracker);
  if (bcf_hdr_write(fout, header) != 0)
    throw VcfWriteException("Failed to write vcf header");

  write_sites(fout, header, gtyper, tracker);

  bcf_close(fout);
}

void populate_vcf_hdr(bcf_hdr_t* hdr, gtyper_ptr gtyper,
                      gram::GenotypeParams const& params,
                      SegmentTracker& tracker) {
  for (auto const& segment : tracker.get_segments()) {
    auto contig =
        vcf_meta_info_line("contig", segment.ID, segment.size).to_string();
    bcf_hdr_append(hdr, contig.c_str());
  }

  auto source = vcf_meta_info_line("source", "gramtools").to_string();
  bcf_hdr_append(hdr, source.c_str());

  bcf_hdr_add_sample(hdr, params.sample_id.c_str());

  for (auto const& header :
       gtyper->get_model_specific_headers()) {  // Genotyper-specific headers
    bcf_hdr_append(hdr, header.to_string().c_str());
  }

  header_vec format_headers = common_headers();  // Common headers
  for (auto& header : format_headers)
    bcf_hdr_append(hdr, header.to_string().c_str());
}

/**
 * Picks only the sites that are not nested in any other (= lvl1 Sites)
 */
std::size_t next_valid_idx(std::size_t idx, std::size_t max_size,
                           gram::parental_map const& p_map) {
  auto result = idx;
  while (result < max_size) {
    auto site_ID = index_to_siteID(result);
    if (p_map.find(site_ID) == p_map.end())
      break;
    else
      result++;
  }
  return result;
}

void write_sites(htsFile* fout, bcf_hdr_t* header, gtyper_ptr const& gtyper,
                 SegmentTracker& tracker) {
  auto p_map = gtyper->get_cov_g()->par_map;
  auto genotyped_records = gtyper->get_genotyped_records();
  std::size_t site_idx{0}, max_size{genotyped_records.size()};
  // Set up and write records
  bcf1_t* record = bcf_init();

  while (true) {
    site_idx = next_valid_idx(site_idx, max_size, p_map);
    if (site_idx >= max_size) break;

    bcf_empty(record);
    populate_vcf_site(header, record, genotyped_records[site_idx], tracker);
    if (bcf_write(fout, header, record) != 0)
      throw VcfWriteException("Failed to write vcf record");
    site_idx++;
  }
}

void add_model_specific_entries(bcf_hdr_t* hdr, bcf1_t* record,
                                site_entries const& entries) {
  for (auto const& entry : entries.doubles) {
    std::vector<float> vals;
    for (auto const& val : entry.vals) vals.push_back((float)val);
    bcf_update_format_float(hdr, record, entry.ID.c_str(), vals.data(),
                            vals.size());
  }
}

void populate_vcf_site(bcf_hdr_t* header, bcf1_t* record, gt_site_ptr site,
                       SegmentTracker& tracker) {
  using str_vec = std::vector<std::string>;

  // Set CHROM
  auto numeric_id =
      bcf_hdr_name2id(header, tracker.get_ID(site->get_pos()).c_str());
  record->rid = numeric_id;
  // Set POS. Pass in a 0-based
  record->pos = tracker.get_relative_pos(site->get_pos());

  // Set GT
  auto gtype_info = site->get_all_gtype_info();
  std::vector<int32_t> gtypes = gtype_info.genotype;
  if (site->is_null()) {
    for (auto& idx : gtypes) idx = bcf_gt_missing;
  } else {
    for (auto& idx : gtypes) idx = bcf_gt_unphased(idx);
  }
  bcf_update_genotypes(header, record, gtypes.data(), gtypes.size());

  // Set DP
  str_vec total_cov{std::to_string(gtype_info.total_coverage)};
  bcf_update_format_string(header, record, "DP",
                           reinterpret_cast<char const**>(total_cov.data()), 1);

  // Set COV
  auto covs = gtype_info.allele_covs;
  if (covs.size() > 0) {
    std::vector<float> depths;
    for (auto const& cov : covs) depths.push_back((float)cov);
    bcf_update_format_float(header, record, "COV", depths.data(),
                            depths.size());
  }

  // Set FT
  str_vec filters = gtype_info.filters;
  if (gtype_info.filters.empty())
    filters.insert(filters.end(), "PASS");
  else {
    auto const max_size = gtype_info.filters.size();
    for (std::size_t i{0}; i < max_size; i++) {
      if (i < max_size - 1) filters.at(i).push_back(',');
    }
  }
  bcf_update_format_string(header, record, "FT",
                           reinterpret_cast<char const**>(filters.data()), 1);

  std::string als;
  for (auto const& al : gtype_info.alleles) {
    if (als.size() != 0) als += ",";
    als += al.sequence;
  }
  bcf_update_alleles_str(header, record, als.c_str());

  add_model_specific_entries(header, record,
                             site->get_model_specific_entries());
}
