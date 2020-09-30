#include "genotype/parameters.hpp"

#include <omp.h>

#include <iostream>

using namespace gram;
using namespace gram::commands::genotype;

struct ploidy_argument {
  Ploidy ploidy;

 public:
  ploidy_argument() = default;
  ploidy_argument(const std::string& in) {
    if (in == "haploid")
      ploidy = Ploidy::Haploid;
    else if (in == "diploid")
      ploidy = Ploidy::Diploid;
    else
      throw std::invalid_argument("Invalid/unsupported ploidy");
  }

  Ploidy get() { return ploidy; }
};

/**
 * Validate ploidy argument
 * Credit:
 * https://www.boost.org/doc/libs/1_72_0/doc/html/program_options/howto.html#id-1.3.32.6.7
 */
void validate(boost::any& v, const std::vector<std::string>& values,
              ploidy_argument* target_type, int) {
  using namespace boost::program_options;

  // Make sure no previous assignment was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  std::string const& s = validators::get_single_string(values);

  v = boost::any(ploidy_argument(s));
}

GenotypeParams commands::genotype::parse_parameters(
    po::variables_map& vm, const po::parsed_options& parsed) {
  GenotypeParams parameters = {};
  std::vector<std::string> reads_fpaths;
  std::string run_dirpath;
  ploidy_argument ploidy;
  Seed::value_type seed;

  po::options_description genotype_description("genotype options");
  genotype_description.add_options()(
      "gram_dir", po::value<std::string>(&parameters.gram_dirpath)->required(),
      "gramtools directory")("reads",
                             po::value<std::vector<std::string>>(&reads_fpaths)
                                 ->multitoken()
                                 ->required(),
                             "file containing reads (FASTA or FASTQ)")(
      "sample_id", po::value<std::string>(&parameters.sample_id)->required())(
      "ploidy", po::value<ploidy_argument>(&ploidy)->required(),
      "expected ploidy of the sample. Choices: {haploid, diploid}")(
      "kmer_size", po::value<uint32_t>(&parameters.kmers_size)->required(),
      "kmer size that got used in build step")(
      "genotype_dir", po::value<std::string>(&run_dirpath)->required(),
      "output directory")("max_threads",
                          po::value<uint32_t>()->default_value(1),
                          "maximum number of threads used")(
      "seed", po::value<SeedSize>(&seed),
      "seed for pseudo-random selection of multi-mapping reads. "
      "a random seed is generated if this option is not used.");

  std::vector<std::string> opts =
      po::collect_unrecognized(parsed.options, po::include_positional);
  if (opts.size() > 0)
    opts.erase(opts.begin());  // Takes out the command itself
  try {
    po::store(po::command_line_parser(opts).options(genotype_description).run(),
              vm);
    po::notify(vm);
  } catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    std::cout << genotype_description << std::endl;
    exit(1);
  }

  fill_common_parameters(parameters, parameters.gram_dirpath);
  for (auto& elem : reads_fpaths) elem = fs::absolute(fs::path(elem)).string();
  parameters.reads_fpaths = reads_fpaths;

  parameters.ploidy = ploidy.get();

  std::string cov_dirpath = mkdir(run_dirpath, "coverage");
  std::string geno_dirpath = mkdir(run_dirpath, "genotype");
  parameters.read_stats_fpath = full_path(run_dirpath, "read_stats.json");
  parameters.debug_fpath =
      full_path(run_dirpath, "site_gtyping_debug_info.txt");

  parameters.allele_sum_coverage_fpath =
      full_path(cov_dirpath, "allele_sum_coverage");
  parameters.allele_base_coverage_fpath =
      full_path(cov_dirpath, "allele_base_coverage.json");
  parameters.grouped_allele_counts_fpath =
      full_path(cov_dirpath, "grouped_allele_counts_coverage.json");

  parameters.genotyped_json_fpath = full_path(geno_dirpath, "genotyped.json");
  parameters.genotyped_vcf_fpath = full_path(geno_dirpath, "genotyped.vcf.gz");
  parameters.personalised_ref_fpath =
      full_path(geno_dirpath, "personalised_reference.fasta");

  parameters.maximum_threads = vm["max_threads"].as<uint32_t>();
  omp_set_num_threads(parameters.maximum_threads);

  if (vm.count("seed")) parameters.seed = seed;
  return parameters;
}