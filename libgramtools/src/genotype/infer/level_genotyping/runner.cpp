#include <cmath>
#include <random>

#include "GCP/GCP.h"

#include "genotype/infer/level_genotyping/runner.hpp"
#include "genotype/infer/allele_extracter.hpp"
#include "prg/coverage_graph.hpp"
#include "prg/make_data_structures.hpp"
#include "genotype/infer/output_specs/fields.hpp"

using namespace gram::genotype::output_spec;

header_vec LevelGenotyper::get_model_specific_headers(){
    header_vec result{
        vcf_meta_info_line{"Model", "LevelGenotyping"},
        vcf_meta_info_line{"FORMAT", "GT_CONF",
                           "Genotype confidence as "
                 "likelihood ratio of called and next most likely genotype.",
                 "1", "Float"},
        vcf_meta_info_line{"FORMAT", "GT_CONF_PERCENTILE",
                    "Percent of calls expected to have lower GT_CONF",
                    "1", "Float"}
    };
    return result;
}

likelihood_related_stats
LevelGenotyper::make_l_stats(double mean_cov_depth, double mean_pb_error){
    PoissonLogPmf poisson_prob{params{mean_cov_depth}};

    // store natural log of pb error also because of its use in likelihood formulae
    return likelihood_related_stats{
        mean_cov_depth,
        mean_pb_error, log(mean_pb_error),
        log(1 - exp(-mean_cov_depth)),
        log(1 - exp(-mean_cov_depth / 2)),
        find_minimum_non_error_cov(mean_pb_error, &poisson_prob),
        poisson_prob,
        PoissonLogPmf{params{mean_cov_depth / 2}}
    };
}

CovCount LevelGenotyper::find_minimum_non_error_cov(double mean_pb_error, poisson_pmf_ptr poisson_prob) {
    CovCount min_count{1};
    while ( (*poisson_prob)(params{static_cast<double>(min_count)}) <= log(pow(mean_pb_error, min_count)) )
        ++min_count;
    return min_count;
}

class ModelDataProducer : public GCP::Model<ModelData>{
private:
    likelihood_related_stats const* l_stats;
    Ploidy ploidy;
public:
    ModelDataProducer(likelihood_related_stats const* l_stats, Ploidy ploidy)
    : l_stats(l_stats), ploidy(ploidy) {};

    ModelData produce_data() override {
        std::binomial_distribution<CovCount> b(l_stats->mean_cov_depth, l_stats->mean_pb_error);
        std::poisson_distribution<CovCount> p(l_stats->mean_cov_depth);
        CovCount const correct_cov = p(random_number_generator);
        CovCount const incorrect_cov = b(random_number_generator);
        allele_vector alleles{
            Allele{"A", {correct_cov}, 0},
            Allele{"A", {incorrect_cov}, 1},
        };
        GroupedAlleleCounts gped_counts{
                {{0}, correct_cov},
                {{1}, incorrect_cov},
        };
        return ModelData(alleles, gped_counts, ploidy, l_stats);
    }
};

/**
 * Draws empirical confidences from the genotyped sites
 * and complements them with simulations if there are not enough.
 */
void LevelGenotyper::compute_confidence_percentiles() {
    constexpr uint16_t distrib_size{10000};
    std::vector<double> confidences(distrib_size);
    auto insertion_point = confidences.begin();

    // Case: draw all needed confidences at random from sites
    if (genotyped_records.size() > distrib_size){
        std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_int_distribution<> distrib(0, genotyped_records.size() - 1);
        while (insertion_point != confidences.end()) {
            auto selected_entry = distrib(generator);
            auto selected_site =
                    std::dynamic_pointer_cast<LevelGenotypedSite>(genotyped_records.at(selected_entry));
            *insertion_point++ = selected_site->get_gt_conf();
        }
    }
    // Case: draw all empirical confidences + simulate up to required number
    else {
        for (auto const& site : genotyped_records)
            *insertion_point++ = std::dynamic_pointer_cast<LevelGenotypedSite>(site)->get_gt_conf();
        auto num_simulations = std::distance(insertion_point, confidences.end());
        // Perform simulations
        GCP::Model<ModelData>* data_producer = new ModelDataProducer(&l_stats, ploidy);
        GCP::Simulator<ModelData, LevelGenotyperModel> simulator(data_producer);
        auto simu = simulator.simulate(num_simulations);
        std::copy(simu.begin(), simu.end(), insertion_point);
    }
    // Add percentiles
    GCP::Percentiler percentiler(confidences);
    for (auto const& site : genotyped_records){
        auto conf = std::dynamic_pointer_cast<LevelGenotypedSite>(site)->get_gt_conf();
        std::dynamic_pointer_cast<LevelGenotypedSite>(site)->set_gt_conf_percentile(
                percentiler.get_confidence_percentile(conf)
                );
    }
}


LevelGenotyper::LevelGenotyper(coverage_Graph const &cov_graph, SitesGroupedAlleleCounts const &gped_covs,
                               ReadStats const &read_stats, Ploidy const ploidy, bool get_gcp) :
                               ploidy(ploidy){
    this->cov_graph = &cov_graph;
    this->gped_covs = &gped_covs;
    child_m = build_child_map(cov_graph.par_map); // Required for site invalidation
    genotyped_records.resize(cov_graph.bubble_map.size()); // Pre-allocate one slot for each bubble in the PRG

    auto mean_cov_depth = read_stats.get_mean_cov_depth();
    auto mean_pb_error = read_stats.get_mean_pb_error();
    l_stats = std::move(make_l_stats(mean_cov_depth, mean_pb_error));

    // Genotype each bubble in the PRG, in most nested to less nested order.
    for (auto const& bubble_pair : cov_graph.bubble_map){
        auto site_ID = bubble_pair.first->get_site_ID();
        auto site_index = siteID_to_index(site_ID);

        auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
        auto extracted_alleles = extracter.get_alleles();
        auto& gped_covs_for_site = gped_covs.at(site_index);

        ModelData data(extracted_alleles, gped_covs_for_site,
                       ploidy, &l_stats, ! extracter.ref_allele_got_made_naturally());
        auto genotyped = LevelGenotyperModel(data);
        auto genotyped_site = genotyped.get_site();
        genotyped_site->set_pos(bubble_pair.first->get_pos());
        genotyped_records.at(site_index) = genotyped_site;
        // Line below is so that when allele extraction occurs and jumps through a previously
        // genotyped site, it knows where in the graph to resume from.
        genotyped_records.at(site_index)->set_site_end_node(bubble_pair.second);

        run_invalidation_process(genotyped_site, site_ID);
    }
    if (get_gcp) compute_confidence_percentiles();
}

