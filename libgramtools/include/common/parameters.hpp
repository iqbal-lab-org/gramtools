/** @file
 * Defines gramtools back-end commands and prg-related filepaths.
 */

#ifndef GRAMTOOLS_PARAMETERS_HPP
#define GRAMTOOLS_PARAMETERS_HPP

#include <string>
#include <vector>

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <filesystem>

namespace po = boost::program_options;
namespace fs = std::filesystem;

namespace gram {

    static std::string mkdir(std::string const& parent_dirpath, std::string const& child_dirpath){
        fs::path dir1(parent_dirpath), dir2(child_dirpath);
        assert(fs::exists(dir1));
        fs::path full_path = dir1 / dir2;
        if (! fs::exists(full_path)) fs::create_directory(full_path);
        return full_path.string();
    }

    // The number of bytes to use for each integer going on disk representing PRG string `Marker`s
    constexpr uint8_t num_bytes_per_integer{4};

    /**
     * PRG file path parameters.
     * Used for either:
     * * serialising all the necessary information for vBWT mapping to a given prg after `build` process.
     * * loading such information for `quasimap`ping reads to the prg.
     * @see gram::commands::build::parse_parameters()
     * @see gram::command::genotype::parse_parameters()
     */
    class CommonParameters {
    public:
        std::string gram_dirpath;
        std::string built_vcf;
        std::string encoded_prg_fpath;
        std::string prg_coords_fpath;
        std::string fm_index_fpath;
        std::string cov_graph_fpath;
        std::string sites_mask_fpath;
        std::string allele_mask_fpath;

        // kmer index file paths
        std::string kmer_index_fpath;
        std::string kmers_fpath;
        std::string kmers_stats_fpath;
        std::string sa_intervals_fpath;
        std::string paths_fpath;

        uint32_t kmers_size;
        uint32_t maximum_threads;
    };

    static std::string full_path(const std::string &base_dirpath,
                                const std::string &file_name) {
        fs::path dir(base_dirpath);
        fs::path file(file_name);
        fs::path full_path = dir / file;
        return full_path.string();
    }

    static void fill_common_parameters(CommonParameters& parameters, std::string const gram_dirpath){
        parameters.gram_dirpath = gram_dirpath;
        parameters.built_vcf = full_path(gram_dirpath, "build.vcf"); // May not exist, if --prg used
        parameters.encoded_prg_fpath = full_path(gram_dirpath, "prg");
        parameters.prg_coords_fpath = full_path(gram_dirpath, "prg_coords.tsv");
        parameters.fm_index_fpath = full_path(gram_dirpath, "fm_index");
        parameters.cov_graph_fpath = full_path(gram_dirpath, "cov_graph");
        parameters.sites_mask_fpath = full_path(gram_dirpath, "variant_site_mask");
        parameters.allele_mask_fpath = full_path(gram_dirpath, "allele_mask");

        parameters.kmer_index_fpath = full_path(gram_dirpath, "kmer_index");
        parameters.kmers_fpath = full_path(gram_dirpath, "kmers");
        parameters.kmers_stats_fpath = full_path(gram_dirpath, "kmers_stats");
        parameters.sa_intervals_fpath = full_path(gram_dirpath, "sa_intervals");
        parameters.paths_fpath = full_path(gram_dirpath, "paths");
    }
}

#endif //GRAMTOOLS_PARAMETERS_HPP
