#include <build/parameters.hpp>
#include "build/build.hpp"


using namespace gram;


void commands::build::run(BuildParams const &parameters) {
    std::cout << "Executing build command" << std::endl;
    auto timer = TimerReport();

    PRG_Info prg_info;

    std::cout << "Loading integer encoded PRG" << std::endl;
    timer.start("Encoded PRG");
    PRG_String ps{parameters.encoded_prg_fpath};
    timer.stop();
    std::cout << "Number of characters in integer encoded linear PRG: "
              << ps.size()
              << std::endl;

    prg_info.last_allele_positions = ps.get_end_positions();

    std::cout << "Generating coverage graph" << std::endl;
    timer.start("Generate Coverage Graph");
    //Need move, not copy assignment, else destructor can affect assigned-to object. Move is compiler default here anyway.
    prg_info.coverage_graph = std::move(generate_cov_graph(parameters, ps));
    timer.stop();

    prg_info.num_variant_sites = prg_info.coverage_graph.bubble_map.size();
    std::cout << "Number of variant sites: " << prg_info.num_variant_sites << std::endl;
    if (prg_info.num_variant_sites == 0) { // No personalised reference to infer; exit
        std::cout << "No variant sites found.\nExiting 1" << std::endl;
        std::exit(1);
    }


    std::cout << "Generating FM-Index" << std::endl;
    timer.start("Generate FM-Index");
    prg_info.fm_index = generate_fm_index(parameters);
    timer.stop();

    std::cout << "Generating PRG masks" << std::endl;
    timer.start("Generating PRG masks");

    prg_info.bwt_markers_mask = generate_bwt_markers_mask(prg_info.fm_index);

    prg_info.dna_bwt_masks = generate_bwt_masks(prg_info.fm_index, parameters);
    prg_info.rank_bwt_a = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_a);
    prg_info.rank_bwt_c = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_c);
    prg_info.rank_bwt_g = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_g);
    prg_info.rank_bwt_t = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_t);
    timer.stop();

    std::cout << "Building kmer index"
              << " (kmer size: " << parameters.kmers_size << ")" << std::endl;
    timer.start("Building kmer index");
    auto kmer_index = kmer_index::build(parameters, prg_info);
    kmer_index::dump(kmer_index, parameters);
    timer.stop();

    timer.report();
}