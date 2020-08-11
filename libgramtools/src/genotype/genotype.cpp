#include "genotype/genotype.hpp"

#include "build/kmer_index/load.hpp"
#include "common/timer_report.hpp"
#include "genotype/infer/level_genotyping/runner.hpp"
#include "genotype/infer/output_specs/make_json.hpp"
#include "genotype/infer/output_specs/make_vcf.hpp"
#include "genotype/infer/output_specs/segment_tracker.hpp"
#include "genotype/infer/personalised_reference.hpp"
#include "genotype/quasimap/quasimap.hpp"

using namespace gram;
using namespace gram::genotype;

namespace gram::genotype {
void write_deduped_p_refs(Fastas const& p_refs, std::string const& fpath) {
  unique_Fastas deduped_p_refs{p_refs.begin(), p_refs.end()};
  std::ofstream pers_ref_fhandle(fpath);
  for (auto& p_ref : deduped_p_refs) pers_ref_fhandle << p_ref << std::endl;
  pers_ref_fhandle.close();
}
}  // namespace gram::genotype

void gram::commands::genotype::run(GenotypeParams const& parameters,
                                   bool const& debug) {
  auto timer = TimerReport();
  /**
   * Quasimap
   */
  std::cout << "Executing genotype command" << std::endl;

  ReadStats readstats;
  std::string first_reads_fpath = parameters.reads_fpaths[0];
  readstats.compute_base_error_rate(first_reads_fpath);

  timer.start("Load data");
  std::cout << "Loading PRG data" << std::endl;
  const auto prg_info = load_prg_info(parameters);
  std::cout << "Loading kmer index data" << std::endl;
  const auto kmer_index = kmer_index::load(parameters);
  timer.stop();

  std::cout << "Running quasimap" << std::endl;
  timer.start("Quasimap");
  auto quasimap_stats =
      quasimap_reads(parameters, kmer_index, prg_info, readstats);

  // Commit the read stats into quasimap output dir.
  std::cout << "Writing read stats to " << parameters.read_stats_fpath
            << std::endl;
  readstats.serialise(parameters.read_stats_fpath);

  std::cout << std::endl;
  std::cout
      << "The following counts include generated reverse complement reads."
      << std::endl;
  std::cout << "Count all reads: " << quasimap_stats.all_reads_count
            << std::endl;
  std::cout << "Count skipped reads: " << quasimap_stats.skipped_reads_count
            << std::endl;
  std::cout << "Count mapped reads: " << quasimap_stats.mapped_reads_count
            << std::endl;
  timer.stop();

  /**
   * Infer
   */
  std::cout << "====================" << std::endl
            << "Running genotyping" << std::endl;
  timer.start("Genotyping");

  // Set up debug info file if needed
  std::string debug_file;
  if (debug) {
    debug_file = parameters.debug_fpath;
    std::cout << "Logging debug genotyping stats to " << debug_file
              << std::endl;
  }

  std::cout << "Running genotyping model" << std::endl;
  LevelGenotyper genotyper{prg_info.coverage_graph,
                           quasimap_stats.coverage.grouped_allele_counts,
                           readstats,
                           parameters.ploidy,
                           true,
                           debug_file};

  std::ifstream coords_file(parameters.prg_coords_fpath);
  SegmentTracker tracker(coords_file);
  coords_file.close();

  std::cout << "Producing json vcf" << std::endl;
  std::ofstream geno_json_fhandle(parameters.genotyped_json_fpath);
  auto gtyper = std::make_shared<LevelGenotyper>(genotyper);
  auto sample_json = make_json_prg(gtyper, tracker);
  sample_json->set_sample_info(parameters.sample_id,
                               "made by gramtools genotype");
  geno_json_fhandle << sample_json->get_prg() << std::endl;
  geno_json_fhandle.close();

  std::cout << "Producing personalised reference" << std::endl;
  auto sites = genotyper.get_genotyped_records();
  tracker.reset();
  auto p_refs =
      get_personalised_ref(prg_info.coverage_graph.root, sites, tracker);
  std::string desc = parameters.sample_id +
                     " personalised reference made by gramtools genotype";
  add_description(p_refs, desc);

  write_deduped_p_refs(p_refs, parameters.personalised_ref_fpath);

  std::cout << "Producing vcf" << std::endl;
  tracker.reset();
  write_vcf(parameters, gtyper, tracker);

  timer.stop();
  timer.report();
}
