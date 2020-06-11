#include <cmath>
#include <random>

#include "GCP/GCP.h"

#include "genotype/infer/level_genotyping/runner.hpp"
#include "genotype/infer/allele_extracter.hpp"
#include "prg/coverage_graph.hpp"
#include "prg/make_data_structures.hpp"
#include "genotype/infer/output_specs/fields.hpp"
#include "genotype/infer/level_genotyping/model.hpp"

using namespace gram::genotype::output_spec;

void add_percentiles(gt_sites const& input_sites, std::vector<double> const& confidences){
    GCP::Percentiler percentiler(confidences);
    for (auto const& site : input_sites){
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

    auto mean_cov = read_stats.get_mean_cov();
    auto var_cov = read_stats.get_var_cov();
    auto mean_pb_error = read_stats.get_mean_pb_error();
    l_stats = make_l_stats(mean_cov, var_cov, mean_pb_error);

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

        // Line below is so that when allele extraction occurs and jumps through a previously
        // genotyped site, it knows where in the graph to resume from.
        genotyped_site->set_site_end_node(bubble_pair.second);

        genotyped_records.at(site_index) = genotyped_site;

        auto downcasted = std::dynamic_pointer_cast<LevelGenotypedSite>(genotyped_site);
        run_invalidation_process(downcasted, site_ID);
        if (genotyped_site->has_filter("AMBIG")) propagate_filter("AMBIG", site_ID);
    }
    if (get_gcp) {
        auto confidences = get_gtconf_distrib(genotyped_records, l_stats, ploidy);
        add_percentiles(genotyped_records, confidences);
    }
}

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

void LevelGenotyper::propagate_filter(std::string const& name, Marker const& parent_site_ID){
    std::vector<Marker> to_process{parent_site_ID};
    Marker cur_ID;
    while (! to_process.empty()){
        cur_ID = to_process.back();
        to_process.pop_back();
        if (child_m.find(cur_ID) == child_m.end()) continue;
        for (auto const& pair: child_m.at(cur_ID)){
            for (auto const& child_marker : pair.second){
                auto referent_site = genotyped_records.at(siteID_to_index(child_marker));
                if (! referent_site->has_filter(name)){
                    referent_site->set_filter(name);
                    to_process.emplace_back(child_marker);
                }
            }
        }
    }
}

void LevelGenotyper::run_invalidation_process(lvlgt_site_ptr const& genotyped_site, Marker const& site_ID) {
    // Invalidation process attempted only if this site contains 1+ site
    if (child_m.find(site_ID) != child_m.end()) {
        auto candidate_haplogroups =
                genotyped_site->get_nonGenotyped_haplogroups();
        auto haplogroups_with_sites = get_haplogroups_with_sites(site_ID, candidate_haplogroups);
        invalidate_if_needed(site_ID, haplogroups_with_sites);
    }
}

AlleleIds LevelGenotyper::get_haplogroups_with_sites(Marker const& site_ID, AlleleIds candidate_haplogroups) const{
    AlleleIds result{};
    if (child_m.find(site_ID) == child_m.end()) return result;
    auto child_entry = child_m.at(site_ID);
    for (auto const& candidate : candidate_haplogroups){
        if (child_entry.find(candidate) != child_entry.end()) result.push_back(candidate);
    }
    return result;
}

void LevelGenotyper::invalidate_if_needed(Marker const& parent_site_ID, AlleleIds haplogroups){
    if (haplogroups.empty()) return;

    std::vector<VariantLocus> to_process;
    for (auto const& haplogroup : haplogroups) to_process.emplace_back(parent_site_ID, haplogroup);

    VariantLocus cur_locus;
    while (! to_process.empty()){
        cur_locus = to_process.back();
        to_process.pop_back();

        // We know the site/haplogroup combination bears 1+ site so use .at() .
        auto sites_on_haplogroup = child_m.at(cur_locus.first).at(cur_locus.second);
        for (auto const& child_marker : sites_on_haplogroup){
            auto referent_site = genotyped_records.at(siteID_to_index(child_marker));
            if (referent_site->is_null()) continue;
            referent_site->make_null();

            auto downcasted = std::dynamic_pointer_cast<LevelGenotypedSite>(referent_site);
            auto all_haplos = downcasted->get_all_haplogroups();
            auto haplos_with_sites = get_haplogroups_with_sites(child_marker, all_haplos);
            for (auto const& h : haplos_with_sites) to_process.emplace_back(child_marker,h);
        }
    }
}


likelihood_related_stats
LevelGenotyper::make_l_stats(double const mean_cov, double const var_cov, double const mean_pb_error) {
    pmf_ptr pmf, pmf_half_depth;
    DataParams data_params{mean_cov, mean_pb_error};
    double prob_no_zero{0}, prob_no_zero_half_depth{0};
    if (var_cov > mean_cov){
        // Case: use negative binomial.
        // Infer its parameters from cov depth mean and variance
        double num_successes = pow(mean_cov, 2) / (var_cov - mean_cov);
        double success_prob = num_successes / (mean_cov + num_successes);
        pmf = std::make_shared<NegBinomLogPmf>(params{num_successes, success_prob});
        prob_no_zero = log(1 - pow(success_prob, num_successes));
        data_params.num_successes = num_successes;
        data_params.success_prob = success_prob;

        num_successes = pow(var_cov, 2) / (var_cov - mean_cov / 2);
        success_prob = num_successes / (mean_cov / 2 + num_successes);
        pmf_half_depth = std::make_shared<NegBinomLogPmf>(params{num_successes, success_prob});
        prob_no_zero_half_depth = log(1 - pow(success_prob, num_successes));
    }
    else {
        pmf = std::make_shared<PoissonLogPmf>(params{mean_cov});
        prob_no_zero = log(1 - exp(mean_cov * -1));

        pmf_half_depth = std::make_shared<PoissonLogPmf>(params{mean_cov / 2});
        prob_no_zero_half_depth = log(1 - exp(mean_cov * -0.5));
    }

    // store natural log of pb error also because of its use in likelihood formulae
    return likelihood_related_stats{
            data_params,
            log(mean_pb_error),
            pmf->operator()(params{0}),
            pmf_half_depth->operator()(params{0}),
            prob_no_zero,
            prob_no_zero_half_depth,
            find_minimum_non_error_cov(mean_pb_error, pmf),
            pmf,
            pmf_half_depth
    };
}

CovCount LevelGenotyper::find_minimum_non_error_cov(double mean_pb_error, pmf_ptr pmf) {
    double min_count{1};
    while ((*pmf)(params{min_count}) <= min_count * log(mean_pb_error))
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
        CovCount correct_cov;
        if (std::dynamic_pointer_cast<PoissonLogPmf>(l_stats->pmf_full_depth)){
            std::poisson_distribution<CovCount> dpois(l_stats->data_params.mean_cov);
            correct_cov = dpois(random_number_generator);
        }
        else{
            std::negative_binomial_distribution<CovCount> dnbinom(
                    l_stats->data_params.num_successes,
                    l_stats->data_params.success_prob);
            correct_cov = dnbinom(random_number_generator);
        }

        std::binomial_distribution<CovCount> b(
                l_stats->data_params.mean_cov,
                l_stats->data_params.mean_pb_error);
        CovCount const incorrect_cov = b(random_number_generator);
        
        allele_vector alleles{
            Allele{"C", {correct_cov}, 0},
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
std::vector<double> LevelGenotyper::get_gtconf_distrib(gt_sites const& input_sites,
                                        likelihood_related_stats const& input_lstats,
                                        Ploidy const& input_ploidy) {
    constexpr uint16_t distrib_size{CONF_DISTRIB_SIZE};
    std::vector<double> confidences(distrib_size);
    auto insertion_point = confidences.begin();

    // Case: draw all needed confidences at random from sites
    if (input_sites.size() > distrib_size){
        std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_int_distribution<> distrib(0, input_sites.size() - 1);
        while (insertion_point != confidences.end()) {
            auto selected_entry = distrib(generator);
            auto selected_site =
                    std::dynamic_pointer_cast<LevelGenotypedSite>(input_sites.at(selected_entry));
            *insertion_point++ = selected_site->get_gt_conf();
        }
    }
    // Case: draw all empirical confidences + simulate up to required number
    else {
        for (auto const& site : input_sites)
            *insertion_point++ = std::dynamic_pointer_cast<LevelGenotypedSite>(site)->get_gt_conf();
        auto num_simulations = std::distance(insertion_point, confidences.end());
        // Perform simulations
        GCP::Model<ModelData>* data_producer = new ModelDataProducer(&input_lstats, input_ploidy);
        GCP::Simulator<ModelData, LevelGenotyperModel> simulator(data_producer);
        auto simu = simulator.simulate(num_simulations);
        std::copy(simu.begin(), simu.end(), insertion_point);
    }
    std::sort(confidences.begin(), confidences.end());
    return confidences;
}


