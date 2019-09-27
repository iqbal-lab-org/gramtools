#include "common/parameters.hpp"
#include "common/timer_report.hpp"

#include "prg/prg.hpp"
#include "prg/load_PRG_string.hpp"
#include "prg/masks.hpp"

#include "kmer_index/build.hpp"
#include "kmer_index/dump.hpp"

#include "build/build.hpp"


using namespace gram;


void commands::build::run(const Parameters &parameters) {
    std::cout << "Executing build command" << std::endl;
    auto timer = TimerReport();

    PRG_Info prg_info;

    std::cout << "Generating integer encoded PRG" << std::endl;
    timer.start("Encoded PRG");
    PRG_String ps{parameters.encoded_prg_fpath};
    prg_info.encoded_prg = ps.get_PRG_string();
    timer.stop();
    std::cout << "Number of characters in integer encoded linear PRG: "
              << prg_info.encoded_prg.size()
              << std::endl;

    prg_info.last_allele_positions = ps.get_end_positions();

    prg_info.max_alphabet_num = get_max_alphabet_num(prg_info.encoded_prg);
    std::cout << "Maximum alphabet character: " << prg_info.max_alphabet_num << std::endl;
    if (prg_info.max_alphabet_num <= 4) { // No personalised reference to infer; exit
        std::cout << "No variant sites found.\nExiting 1" << std::endl;
        std::exit(1);
    }

    std::cout << "Generating FM-Index" << std::endl;
    timer.start("Generate FM-Index");
    prg_info.fm_index = generate_fm_index(parameters);
    timer.stop();

    std::cout << "Generating PRG masks" << std::endl;
    timer.start("Generating PRG masks");
    prg_info.dna_bwt_masks = generate_bwt_masks(prg_info.fm_index, parameters);

    prg_info.sites_mask = generate_sites_mask(prg_info.encoded_prg);
    sdsl::store_to_file(prg_info.sites_mask, parameters.sites_mask_fpath);

    prg_info.allele_mask = generate_allele_mask(prg_info.encoded_prg);
    sdsl::store_to_file(prg_info.allele_mask, parameters.allele_mask_fpath);

    prg_info.prg_markers_mask = generate_prg_markers_mask(prg_info.encoded_prg);
    prg_info.prg_markers_rank = sdsl::rank_support_v<1>(&prg_info.prg_markers_mask);
    prg_info.prg_markers_select = sdsl::select_support_mcl<1>(&prg_info.prg_markers_mask);

    prg_info.markers_mask_count_set_bits =
            prg_info.prg_markers_rank(prg_info.prg_markers_mask.size());

    prg_info.bwt_markers_mask = generate_bwt_markers_mask(prg_info.fm_index);

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