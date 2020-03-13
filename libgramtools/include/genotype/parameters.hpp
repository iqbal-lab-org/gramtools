/**
 * @file
 * Command-line argument processing for `quasimap` command.
 */


#ifndef GRAMTOOLS_QUASIMAP_PARAMETERS_HPP
#define GRAMTOOLS_QUASIMAP_PARAMETERS_HPP

#include "common/parameters.hpp"

namespace gram{
    enum class Ploidy{Haploid, Diploid};

    class GenotypeParams : public CommonParameters {
    public:
        std::vector<std::string> reads_fpaths;

        std::string allele_sum_coverage_fpath;
        std::string allele_base_coverage_fpath;
        std::string grouped_allele_counts_fpath;
        std::string read_stats_fpath;

        Ploidy ploidy;
        std::string sample_id;
        std::string genotyped_json_fpath;
        std::string genotyped_vcf_fpath;
        std::string personalised_ref_fpath;

        uint32_t seed;
    };

    namespace commands::genotype {
    /**
     * Parse command line parameters.
     * A directory containing the information necessary for vBWT mapping to the prg must be passed.
     * This contains files produced in the `build` stage. They are then loaded during `quasimap`ping.
     * @see gram::commands::build::run()
     */
    GenotypeParams parse_parameters(po::variables_map &vm,
                                    const po::parsed_options &parsed);
    }
}

#endif //GRAMTOOLS_QUASIMAP_PARAMETERS_HPP
