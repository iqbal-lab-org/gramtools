#include <iomanip>

#include "prg/coverage_graph.hpp"
#include "prg/make_data_structures.hpp"
#include "genotype/infer/output_specs/make_json.hpp"
#include "genotype/infer/output_specs/segment_tracker.hpp"
#include "genotype/infer/allele_extracter.hpp"
#include "genotype/infer/personalised_reference.hpp"

#include "simulate/simulate.hpp"
#include "simulate/induce_genotypes.hpp"

using namespace gram::genotype;
using namespace gram::simulate;

namespace gram::simulate {

RandomGenotyper::RandomGenotyper(coverage_Graph const &cov_graph) {
        this->cov_graph = &cov_graph;
        child_m = build_child_map(cov_graph.par_map); // Required for site invalidation
        genotyped_records.resize(cov_graph.bubble_map.size()); // Pre-allocate one slot for each bubble in the PRG

        // Genotype each bubble in the PRG, in most nested to less nested order.
        for (auto const &bubble_pair : cov_graph.bubble_map) {
            auto site_ID = bubble_pair.first->get_site_ID();
            auto site_index = siteID_to_index(site_ID);

            auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
            RandomInclusiveInt rand(0);
            auto genotyped_site = make_randomly_genotyped_site(&rand, extracter.get_alleles(),
                    extracter.ref_allele_got_made_naturally());
            genotyped_site->set_pos(bubble_pair.first->get_pos());
            genotyped_site->set_site_end_node(bubble_pair.second);

            genotyped_records.at(site_index) = genotyped_site;

            run_invalidation_process(genotyped_site, site_ID);
        }
    }

    gt_site_ptr make_randomly_genotyped_site(RandomGenerator const *const rand, allele_vector const &alleles,
                                             bool const use_ref_allele) {
        allele_vector picked_alleles{alleles.begin(), alleles.begin() + 1};  // Always pick REF
        uint32_t picked_index;
        if (use_ref_allele) picked_index = rand->generate(0, alleles.size() - 1);
        else picked_index = rand->generate(1, alleles.size() - 1);
        auto chosen_hapg = alleles.at(picked_index).haplogroup;
        allele_coverages covs{1};

        if (picked_index != 0) {
            picked_alleles.push_back(alleles.at(picked_index));
            covs.push_back(1);
            covs.at(0) = 0;
            picked_index = 1;
        }

        auto result = std::make_shared<SimulatedSite>();
        result->populate_site(gtype_information{
           picked_alleles,
           GtypedIndices{static_cast<int>(picked_index)},
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

    std::stringstream coords_file{""};
    SegmentTracker tracker(coords_file);

    uint64_t num_runs{0};
    uint64_t num_sampled{0};
    while (num_runs < parameters.max_num_paths){
       RandomGenotyper gtyper(cov_g);
       auto ptr_gtyper = std::make_shared<RandomGenotyper>(gtyper);
       auto genotyped_records = gtyper.get_genotyped_records();
       auto new_p_ref = get_personalised_ref(cov_g.root, genotyped_records, tracker).at(0);

       if (unique_simu_paths.find(new_p_ref) == unique_simu_paths.end()){
           num_sampled++;
           sample_id = parameters.sample_id + std::to_string(num_sampled);
           new_p_ref.set_ID(sample_id);
           new_p_ref.set_desc("made by gramtools simulate");
           unique_simu_paths.insert(new_p_ref);
           ordered_simu_paths.push_back(new_p_ref);
           auto new_json = make_json_prg(ptr_gtyper, tracker);
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

header_vec RandomGenotyper::get_model_specific_headers(){
   return header_vec{
       vcf_meta_info_line{
           "Model", "Simulated_Path"
       }
   };
}
