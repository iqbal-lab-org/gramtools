#include <htslib/synced_bcf_reader.h>
#include "genotype/infer/output_specs/make_vcf.hpp"
#include "genotype/infer/output_specs/fields.hpp"
#include "prg/coverage_graph.hpp"
#include <filesystem>

namespace fs = std::filesystem;

void write_vcf(gram::GenotypeParams const &params, gtyper_ptr const &gtyper) {
    auto fout = bcf_open(params.genotyped_vcf_fpath.c_str(), "wz"); // Writer
    bcf_srs_t* sr = nullptr; // Reader, if a vcf was used to build the PRG

    // Set up and write header
    bcf_hdr_t* hdr;
    if (fs::exists(params.built_vcf)) {
        sr = bcf_sr_init();
        bcf_sr_add_reader(sr, params.built_vcf.c_str());
        hdr = bcf_sr_get_header(sr, 0);
    }
    else hdr = bcf_hdr_init("w");
    populate_vcf_hdr(hdr, gtyper, params);
    bcf_hdr_write(fout, hdr);

    write_sites(fout, sr, hdr, gtyper);

    if (sr != nullptr){
        if ( sr->errnum ) ("Error: %s\n", bcf_sr_strerror(sr->errnum));
        bcf_sr_destroy(sr);
    }

    bcf_close(fout);
}

void populate_vcf_hdr(bcf_hdr_t *hdr, gtyper_ptr gtyper,
        gram::GenotypeParams const &params) {

    if (fs::exists(params.built_vcf)){
        bcf_hdr_remove(hdr, 2, nullptr); // Remove all FORMAT header entries
        bcf_hdr_remove(hdr, 5, "source");
        bcf_hdr_set_samples(hdr, nullptr, 0); // Remove samples, from reading and writing
    }
    else {
       auto contig = vcf_meta_info_line("contig", "gram_prg", "").to_string();
       bcf_hdr_append(hdr, contig.c_str());
    }

    auto source = vcf_meta_info_line("source", "gramtools").to_string();
    bcf_hdr_append(hdr, source.c_str());

    bcf_hdr_add_sample(hdr, params.sample_id.c_str());

    for (auto const& header : gtyper->get_model_specific_headers()){ // Genotyper-specific headers
        bcf_hdr_append(hdr, header.to_string().c_str());
    }

    header_vec format_headers = vcf_format_headers(); // Common headers
    for (auto& header: format_headers)
        bcf_hdr_append(hdr, header.to_string().c_str());
}

/**
 * Picks only the sites that are not nested in any other (= lvl1 Sites)
 */
std::size_t next_valid_idx(std::size_t idx, std::size_t max_size, gram::parental_map const& p_map){
   auto result = idx;
   while (result < max_size){
      auto site_ID = index_to_siteID(result);
      if (p_map.find(site_ID) == p_map.end()) break;
      else result++;
   }
   return result;
}


void write_sites(htsFile *fout, bcf_srs_t *sr, bcf_hdr_t *hdr, gtyper_ptr const &gtyper) {
    auto p_map = gtyper->get_cov_g()->par_map;
    auto genotyped_records = gtyper->get_genotyped_records();
    std::size_t site_idx{0}, max_size{genotyped_records.size()};
    // Set up and write records
    bcf1_t* record = bcf_init();

    while (true){
        site_idx = next_valid_idx(site_idx, max_size, p_map);
        if (site_idx >= max_size) break;

        if (sr != nullptr){
            if (bcf_sr_next_line(sr) == 0) break;
            else record = bcf_sr_get_line(sr, 0);
        }
        else bcf_empty(record);

        populate_vcf_site(hdr, record, genotyped_records[site_idx]);
        bcf_write(fout, hdr, record);
        site_idx++;
    }
}

void add_model_specific_entries(bcf_hdr_t* hdr, bcf1_t* record, site_entries const& entries){
    for (auto const& entry: entries.doubles){
        std::vector<float> vals;
        for (auto const& val : entry.vals) vals.push_back((float)val);
        bcf_update_format_float(hdr, record, entry.ID.c_str(), vals.data(), vals.size());
    }
}

void populate_vcf_site(bcf_hdr_t* hdr, bcf1_t* record, gt_site_ptr site){
    using str_vec = std::vector<std::string>;

    record->pos = site->get_pos();
    auto gtype_info = site->get_all_gtype_info();
    std::vector<int32_t> gtypes = gtype_info.genotype;
    if (site->is_null()){
        for (auto& idx : gtypes) idx = bcf_gt_missing;
    }
    else{
        for (auto& idx : gtypes) idx = bcf_gt_unphased(idx);
    }
    bcf_update_genotypes(hdr, record, gtypes.data(), gtypes.size());

    str_vec total_cov{std::to_string(gtype_info.total_coverage)};
    bcf_update_format_string(hdr, record, "DP", reinterpret_cast<char const **>(total_cov.data()), 1);

    auto covs = gtype_info.allele_covs;
    if (covs.size() > 0){
        std::vector<float> depths;
        for (auto const& cov : covs) depths.push_back((float)cov);
        bcf_update_format_float(hdr, record, "COV", depths.data(), depths.size());
    }

    std::string als;
    for (auto const& al : gtype_info.alleles) {
        if (als.size() != 0) als += ",";
        als += al.sequence;
    }
    bcf_update_alleles_str(hdr, record, als.c_str());

    add_model_specific_entries(hdr, record, site->get_model_specific_entries());
}
