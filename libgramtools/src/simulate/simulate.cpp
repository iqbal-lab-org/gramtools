#include "simulate/simulate.hpp"
#include "prg/coverage_graph.hpp"
#include "genotype/infer/output_specs/json_prg_spec.hpp"
#include "genotype/infer/output_specs/json_site_spec.hpp"
#include "genotype/infer/allele_extracter.hpp"
#include "genotype/infer/personalised_reference.hpp"
#include <iomanip>

using namespace gram::genotype;
using namespace gram::simulate;

namespace gram::simulate {

    RandomGenotypedSite::RandomGenotypedSite(){
        auto json_site_ptr = std::make_shared<Json_Site>();
        json_site = json_site_ptr;
    }

RandomGenotyper::RandomGenotyper(coverage_Graph const &cov_graph) {
        this->cov_graph = &cov_graph;
        child_m = build_child_map(cov_graph.par_map); // Required for site invalidation
        genotyped_records.resize(cov_graph.bubble_map.size()); // Pre-allocate one slot for each bubble in the PRG

        auto json_ptr = std::make_shared<Simulated_Json>();
        this->json_prg = json_ptr;

        // Genotype each bubble in the PRG, in most nested to less nested order.
        for (auto const &bubble_pair : cov_graph.bubble_map) {
            auto site_ID = bubble_pair.first->get_site_ID();
            auto site_index = siteID_to_index(site_ID);

            auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
            RandomInclusiveInt rand(0);
            auto genotyped_site = make_randomly_genotyped_site(&rand, extracter.get_alleles());

            genotyped_records.at(site_index) = genotyped_site;
            // Line below is so that when allele extraction occurs and jumps through a previously
            // genotyped site, it knows where in the graph to resume from.
            genotyped_records.at(site_index)->set_site_end_node(bubble_pair.second);

            run_invalidation_process(genotyped_site, site_ID);
        }
    }

    gt_site_ptr make_randomly_genotyped_site(RandomGenerator const *const rand, allele_vector const &alleles) {
        allele_vector picked_alleles{alleles.begin(), alleles.begin() + 1};  // Always pick REF
        auto picked_index = rand->generate(0, alleles.size() - 1);
        auto chosen_hapg = alleles.at(picked_index).haplogroup;
        allele_coverages covs{1};

        if (picked_index != 0) {
            picked_alleles.push_back(alleles.at(picked_index));
            covs.push_back(1);
            covs.at(0) = 0;
            picked_index = 1;
        }

        auto result = std::make_shared<RandomGenotypedSite>();
        result->populate_site(gtype_information{
           picked_alleles,
           GtypedIndices{picked_index},
           covs,
           1,
           AlleleIds{chosen_hapg}
        });
        // Note: each picked allele is from a different haplogroup, because we always pick a single genotyped allele
        result->set_num_haplogroups(alleles.size()); // Needed for invalidation process

        return result;
    }
}

void gram::commands::simulate::run(SimulateParams const& parameters){
    PRG_String prg{parameters.encoded_prg_fpath};
    coverage_Graph cov_g{prg};

    std::string desc{"path through prg made by gramtools simulate"};
    unique_Fastas unique_simu_paths;
    Fastas ordered_simu_paths;
    json_prg_ptr simu_json;
    bool first{true};
    std::string sample_id;

    uint64_t num_runs{0};
    uint64_t num_sampled{0};
    while (num_runs < parameters.max_num_paths){
       RandomGenotyper gtyper(cov_g);
       auto genotyped_records = gtyper.get_genotyped_records();
       auto new_p_ref = get_personalised_ref(cov_g.root, genotyped_records).at(0);

       if (unique_simu_paths.find(new_p_ref) == unique_simu_paths.end()){
           num_sampled++;
           sample_id = parameters.sample_id + std::string("_") + std::to_string(num_sampled);
           new_p_ref.set_sample_info(sample_id, desc);
           unique_simu_paths.insert(new_p_ref);
           ordered_simu_paths.push_back(new_p_ref);
           auto new_json = gtyper.get_JSON();
           if (first){
               simu_json = new_json;
               simu_json->set_sample_info(sample_id, desc);
               first = false;
           }
           else {
               new_json->set_sample_info(sample_id, desc);
               simu_json->combine_with(*new_json);
           }
       }
        num_runs++;
    }

    std::cout << "Made " << unique_simu_paths.size() << " simulated paths." << std::endl;

    // Write fasta
    std::ofstream fasta_fhandle(parameters.fasta_out_fpath);
    for (auto& simu_path : ordered_simu_paths) fasta_fhandle << simu_path << std::endl;
    fasta_fhandle.close();

    // Write JSON
    std::ofstream geno_json_fhandle(parameters.json_out_fpath);
    geno_json_fhandle << std::setw(4) << simu_json->get_prg() << std::endl;
    geno_json_fhandle.close();
}
