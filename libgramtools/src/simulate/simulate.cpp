#include <boost/iostreams/filtering_streambuf.hpp>

#include "prg/coverage_graph.hpp"
#include "prg/make_data_structures.hpp"
#include "genotype/infer/output_specs/make_json.hpp"
#include "genotype/infer/output_specs/segment_tracker.hpp"
#include "genotype/infer/allele_extracter.hpp"
#include "genotype/infer/personalised_reference.hpp"
#include "common/file_read.hpp"

#include "simulate/simulate.hpp"
#include "simulate/induce_genotypes.hpp"


using namespace gram::genotype;
using namespace gram::simulate;

namespace gram::simulate {

SimulationGenotyper::SimulationGenotyper(coverage_Graph const &cov_graph) {
        this->cov_graph = &cov_graph;
        child_m = build_child_map(cov_graph.par_map); // Required for site invalidation & json output
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

    lvlgt_site_ptr make_randomly_genotyped_site(RandomGenerator const *const rand, allele_vector const &alleles,
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

SimulationGenotyper::SimulationGenotyper(coverage_Graph const &cov_graph, gt_sites const& input_sites){
    this->cov_graph = &cov_graph;
    child_m = build_child_map(cov_graph.par_map); // Required for site invalidation & json output
    genotyped_records = input_sites;
}

void add_new_json(json_prg_ptr& simu_json, gtyper_ptr const& gtyper, SegmentTracker& tracker,
        bool& first, std::string const& sample_id, std::string const& desc){
    auto new_json = make_json_prg(gtyper, tracker);
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

void simulate_paths(json_prg_ptr& simu_json, coverage_Graph const& cov_graph,
                    SimulateParams const& parameters){
    std::string desc{"path through prg made by gramtools simulate"}, sample_id;
    unique_Fastas unique_simu_paths;
    Fastas ordered_simu_paths;
    bool first{true};

    std::stringstream coords_file{""};
    SegmentTracker tracker(coords_file);

    uint64_t num_runs{0}, num_sampled{0};
    while (num_runs < parameters.max_num_paths){
        num_runs++;
        auto gtyper = std::make_shared<SimulationGenotyper>(cov_graph);
        auto genotyped_records = gtyper->get_genotyped_records();
        auto new_p_ref = get_personalised_ref(cov_graph.root, genotyped_records, tracker).at(0);

        if (unique_simu_paths.find(new_p_ref) != unique_simu_paths.end()) continue;

        num_sampled++;
        sample_id = parameters.sample_id + std::to_string(num_sampled);
        new_p_ref.set_ID(sample_id);
        new_p_ref.set_desc("made by gramtools simulate");
        unique_simu_paths.insert(new_p_ref);
        ordered_simu_paths.push_back(new_p_ref);

        add_new_json(simu_json, gtyper, tracker, first, sample_id, desc);
    }

    std::cout << "Made " << unique_simu_paths.size() << " simulated paths." << std::endl;

    // Write fasta
    std::ofstream fasta_fhandle(parameters.fasta_out_fpath);
    for (auto& simu_path : ordered_simu_paths) fasta_fhandle << simu_path << std::endl;
    fasta_fhandle.close();
}


gt_sites induce_genotypes_one_seq(gt_sites const& template_sites, coverage_Graph const& input_prg,
        std::string const& sequence, std::string const& seq_id){
   gt_sites result;
   result.reserve(template_sites.size());
   for (auto const& template_site : template_sites){
       auto downcasted = std::dynamic_pointer_cast<SimulatedSite>(template_site);
       result.push_back(std::make_shared<SimulatedSite>(*downcasted));
   }

   nt_ptr endpoint;
   endpoint = thread_sequence(input_prg.root, sequence, seq_id, false);

   apply_genotypes(endpoint, result);
   return result;
}

void induce_genotypes_all_seqs(json_prg_ptr& simu_json, coverage_Graph const& input_prg,
        std::string const& fasta_fpath){
    boost::iostreams::filtering_istreambuf in;
    std::ifstream fhandle(fasta_fpath, std::ios::binary);
    if (! fhandle.is_open())
        throw std::ios_base::failure("Could not open: " + fasta_fpath);
    input_fasta(in, fhandle, is_gzipped(fasta_fpath));
    std::istream getter{&in};

    auto template_sites = make_nulled_sites(input_prg);
    bool first{true};
    std::stringstream coords_file{""};
    SegmentTracker tracker(coords_file);
    std::string line, fasta_id, fasta_seq, desc{"induced genotypes made by gramtools simulate"};
    while(std::getline(getter, line)) {
        if (line[0] == '>') {
            if (!fasta_seq.empty()) {
                auto gtyped_sites = induce_genotypes_one_seq(template_sites, input_prg, fasta_seq, fasta_id);
                auto gtyper = std::make_shared<SimulationGenotyper>(input_prg, gtyped_sites);
                add_new_json(simu_json, gtyper, tracker, first, fasta_id, desc);
            }
            fasta_seq.clear();
            fasta_id = line.substr(1, line.find(" ") - 1); // Take first word of header only
        }
        else if (! fasta_id.empty())
                fasta_seq += line;
    }
    if (! fasta_seq.empty() && !fasta_id.empty()) {
        auto gtyped_sites = induce_genotypes_one_seq(template_sites, input_prg, fasta_seq, fasta_id);
        auto gtyper = std::make_shared<SimulationGenotyper>(input_prg, gtyped_sites);
        add_new_json(simu_json, gtyper, tracker, first, fasta_id, desc);
    }
}

void gram::commands::simulate::run(SimulateParams const& parameters){
    PRG_String prg{parameters.encoded_prg_fpath};
    coverage_Graph cov_g{prg};
    json_prg_ptr simu_json;

    if (parameters.input_sequences_fpath.empty())
       simulate_paths(simu_json, cov_g, parameters);
    else induce_genotypes_all_seqs(simu_json, cov_g, parameters.input_sequences_fpath);

    // Write JSON
    std::ofstream geno_json_fhandle(parameters.json_out_fpath);
    geno_json_fhandle << simu_json->get_prg() << std::endl;
    geno_json_fhandle.close();
}

header_vec SimulationGenotyper::get_model_specific_headers(){
   return header_vec{
       vcf_meta_info_line{
           "Model", "Simulated_Path"
       }
   };
}
