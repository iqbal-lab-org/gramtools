#include <htslib/synced_bcf_reader.h>
#include <boost/filesystem.hpp>
#include "genotype/infer/output_specs/make_vcf.hpp"
#include "genotype/infer/output_specs/fields.hpp"

namespace fs = boost::filesystem;

void write_vcf(gram::GenotypeParams const& params, gt_sites const& sites,
               gtyper_ptr const& gtyper){
    auto fout = bcf_open(params.genotyped_vcf_fpath.c_str(), "wz");

    bcf_srs_t* sr = nullptr;
    bcf_hdr_t* hdr;
    if (fs::exists(params.built_vcf)) {
        sr = bcf_sr_init();
        bcf_sr_add_reader(sr, params.built_vcf.c_str());
        hdr = bcf_sr_get_header(sr, 0);
    }
    else hdr = bcf_hdr_init("w");
    populate_vcf_hdr(hdr, gtyper, params);
    bcf_hdr_write(fout, hdr);

    bcf1_t* record;
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
    }

    bcf_hdr_append(hdr,
                   vcf_meta_info_line("source", "gramtools").c_str());
    bcf_hdr_add_sample(hdr, params.sample_id.c_str());

    headers format_headers = vcf_format_headers();
    for (auto& entry: format_headers)
        bcf_hdr_append(hdr, entry.second.c_str());
}
