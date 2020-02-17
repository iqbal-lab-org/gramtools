#include "genotype/genotype.hpp"
#include "genotype/quasimap/quasimap.hpp"
#include "genotype/infer/level_genotyping/runner.hpp"
#include "genotype/infer/personalised_reference.hpp"
#include "genotype/infer/json_spec/prg_spec.hpp"

#include "common/timer_report.hpp"
#include "kmer_index/load.hpp"

using namespace gram;
using namespace gram::genotype;

void gram::commands::genotype::run(GenotypeParams const& parameters){
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
    auto quasimap_stats = quasimap_reads(parameters, kmer_index, prg_info, readstats);

    // Commit the read stats into quasimap output dir.
    std::cout << "Writing read stats to " << parameters.read_stats_fpath << std::endl;
    readstats.serialise(parameters.read_stats_fpath);

    std::cout << std::endl;
    std::cout << "The following counts include generated reverse complement reads."
              << std::endl;
    std::cout << "Count all reads: " << quasimap_stats.all_reads_count << std::endl;
    std::cout << "Count skipped reads: " << quasimap_stats.skipped_reads_count << std::endl;
    std::cout << "Count mapped reads: " << quasimap_stats.mapped_reads_count << std::endl;
    timer.stop();


    /**
     * Infer
     */
    std::cout << "Running genotyping" << std::endl;
    timer.start("Genotyping");
    LevelGenotyper genotyper{prg_info.coverage_graph, quasimap_stats.coverage.grouped_allele_counts,
                             readstats, parameters.ploidy};
    timer.stop();
    std::ofstream geno_json_fhandle(parameters.genotyped_json_fpath);
    auto sample_json = genotyper.get_JSON();
    //sample_json->set_sample_info();
    geno_json_fhandle << std::setw(4) << sample_json << std::endl;
    geno_json_fhandle.close();

    auto sites = genotyper.get_genotyped_records();
    auto p_refs = get_personalised_ref(prg_info.coverage_graph.root, sites);
    std::ofstream pers_ref_fhandle(parameters.personalised_ref_fpath);
    for (auto const& p_ref : p_refs) pers_ref_fhandle << p_ref << std::endl;
    pers_ref_fhandle.close();

    timer.report();
}

