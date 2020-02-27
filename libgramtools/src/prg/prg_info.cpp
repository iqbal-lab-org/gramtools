#include "build/kmer_index/masks.hpp"
#include "prg/prg_info.hpp"


using namespace gram;

PRG_Info gram::load_prg_info(CommonParameters const &parameters) {
    PRG_Info prg_info;

    PRG_String ps{parameters.encoded_prg_fpath};
    prg_info.last_allele_positions = ps.get_end_positions();

    // Load coverage graph
    std::ifstream ifs{parameters.cov_graph_fpath};
    boost::archive::binary_iarchive ia{ifs};
    ia >> prg_info.coverage_graph;
    prg_info.num_variant_sites = prg_info.coverage_graph.bubble_map.size();

    prg_info.fm_index = load_fm_index(parameters);

    prg_info.bwt_markers_mask = generate_bwt_markers_mask(prg_info.fm_index);

    prg_info.dna_bwt_masks = load_dna_bwt_masks(prg_info.fm_index, parameters);
    prg_info.rank_bwt_a = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_a);
    prg_info.rank_bwt_c = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_c);
    prg_info.rank_bwt_g = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_g);
    prg_info.rank_bwt_t = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_t);

    return prg_info;
}

