#include "common/utils.hpp"
#include "genotype/quasimap/quasimap.hpp"
#include "genotype/parameters.hpp"
#include <boost/filesystem.hpp>


using namespace gram;
namespace fs = boost::filesystem;

struct ploidy_argument{
    Ploidy ploidy;
public:
    ploidy_argument() = default;
    ploidy_argument(const std::string& in){
       if (in == "haploid") ploidy = Ploidy::Haploid;
       else if (in == "diploid") ploidy = Ploidy::Diploid;
       else throw std::invalid_argument("Invalid/unsupported ploidy");
    }

    Ploidy get() {return ploidy;}
};


/**
 * Validate ploidy argument
 * Credit: https://www.boost.org/doc/libs/1_72_0/doc/html/program_options/howto.html#id-1.3.32.6.7
 */
void validate(boost::any& v,
              const std::vector<std::string>& values,
              ploidy_argument* target_type, int)
{
    using namespace boost::program_options;

    // Make sure no previous assignment was made.
    validators::check_first_occurrence(v);
    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    std::string const& s = validators::get_single_string(values);

    v = boost::any(ploidy_argument(s));
}

Parameters commands::genotype::parse_parameters(po::variables_map &vm,
                                                const po::parsed_options &parsed) {
    std::string gram_dirpath;
    std::string run_dirpath;
    ploidy_argument ploidy;
    uint32_t kmer_size;
    std::vector<std::string> reads;

    po::options_description genotype_description("genotype options");
    genotype_description.add_options()
                                ("gram_dir", po::value<std::string>(&gram_dirpath)->required(),
                                 "gramtools directory")
                                ("reads", po::value<std::vector<std::string>>(&reads)->multitoken()->required(),
                                 "file containing reads (FASTA or FASTQ)")
                                ("ploidy", po::value<ploidy_argument>(&ploidy)->required(),
                                        "expected ploidy of the sample. Choices: {haploid, diploid}")
                                ("kmer_size", po::value<uint32_t>(&kmer_size)->required(),
                                 "kmer size that got used in build step")
                                ("genotype_dir", po::value<std::string>(&run_dirpath)->required(),
                                 "output directory")
                                ("max_threads", po::value<uint32_t>()->default_value(1),
                                 "maximum number of threads used")
                                ("seed", po::value<uint32_t>()->default_value(0),
                                        "seed for pseudo-random selection of multi-mapping reads. "
                                        "the default of 0 produces a random seed.");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    if (opts.size() > 0) opts.erase(opts.begin()); // Takes out the command itself
    try {
        po::store(po::command_line_parser(opts).options(genotype_description).run(), vm);
        po::notify(vm);
    }
    catch(const std::exception &e){
        std::cout << e.what() << std::endl;
        std::cout << genotype_description << std::endl;
        exit(1);
    }

    Parameters parameters = {};
    parameters.gram_dirpath = gram_dirpath;
    parameters.encoded_prg_fpath = full_path(gram_dirpath, "prg");
    parameters.fm_index_fpath = full_path(gram_dirpath, "fm_index");
    parameters.cov_graph_fpath = full_path(gram_dirpath, "cov_graph");
    parameters.sites_mask_fpath = full_path(gram_dirpath, "variant_site_mask");
    parameters.allele_mask_fpath = full_path(gram_dirpath, "allele_mask");

    parameters.kmer_index_fpath = full_path(gram_dirpath, "kmer_index");
    parameters.kmers_fpath = full_path(gram_dirpath, "kmers");
    parameters.kmers_stats_fpath = full_path(gram_dirpath, "kmers_stats");
    parameters.sa_intervals_fpath = full_path(gram_dirpath, "sa_intervals");
    parameters.paths_fpath = full_path(gram_dirpath, "paths");

    parameters.ploidy = ploidy.get();
    parameters.kmers_size = kmer_size;
    parameters.reads_fpaths = reads;

    std::string cov_dirpath = mkdir(run_dirpath, "coverage");
    std::string geno_dirpath = mkdir(run_dirpath, "genotype");
    parameters.sdsl_memory_log_fpath = full_path(run_dirpath, "sdsl_memory_log");
    parameters.read_stats_fpath = full_path(run_dirpath, "read_stats.json");

    parameters.allele_sum_coverage_fpath = full_path(cov_dirpath, "allele_sum_coverage");
    parameters.allele_base_coverage_fpath = full_path(cov_dirpath, "allele_base_coverage.json");
    parameters.grouped_allele_counts_fpath = full_path(cov_dirpath, "grouped_allele_counts_coverage.json");

    parameters.genotyped_json = full_path(geno_dirpath, "genotyped.json");

    parameters.maximum_threads = vm["max_threads"].as<uint32_t>();
    parameters.seed = vm["seed"].as<uint32_t>();
    return parameters;
}