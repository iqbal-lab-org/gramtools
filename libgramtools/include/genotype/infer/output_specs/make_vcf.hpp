#ifndef MAKE_VCF_HPP
#define MAKE_VCF_HPP

#include <htslib/vcf.h>
#include "genotype/infer/interfaces.hpp"
#include "genotype/parameters.hpp"

using namespace gram::genotype::infer;

void populate_vcf_hdr(bcf_hdr_t *hdr, gtyper_ptr gtyper, gram::GenotypeParams const &params);
void write_vcf(gram::GenotypeParams const& params, gt_sites const& sites,
        gtyper_ptr const& gtyper);

#endif //MAKE_VCF_HPP
