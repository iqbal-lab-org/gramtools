#include "common/utils.hpp"
#include "kmer_index/masks.hpp"
#include "common.hpp"


using namespace gram;


std::string gram::submods::decode(const uint64_t base){
    switch(base){
        case 1: return "A";
        case 2: return "C";
        case 3: return "G";
        case 4: return "T";
        default: return std::to_string(base);
    }
}


PRG_Info gram::submods::generate_prg_info(const marker_vec &prg_raw) {
    BuildParams parameters = {};
   parameters.encoded_prg_fpath = "encoded_prg_file_name";
    parameters.fm_index_fpath = "@fm_index";
    parameters.gram_dirpath = "@gram_dir";


    PRG_String ps{prg_raw};
    auto encoded_prg = ps.get_PRG_string();
    // Write the int vector to disk so that it can be read by sdsl for building fm index
    ps.write(parameters.encoded_prg_fpath, endianness::little);

    PRG_Info prg_info;
    prg_info.encoded_prg = encoded_prg;
    prg_info.fm_index = generate_fm_index(parameters);
    // NB: the move is crucial here, otherwise the initialised cov_Graph's destructor affects the assigned-to cov_Graph
    prg_info.coverage_graph = std::move(coverage_Graph{ps});
    prg_info.last_allele_positions = ps.get_end_positions();
    prg_info.sites_mask = generate_sites_mask(encoded_prg);
    prg_info.allele_mask = generate_allele_mask(encoded_prg);

    prg_info.prg_markers_mask = generate_prg_markers_mask(encoded_prg);
    prg_info.prg_markers_rank = sdsl::rank_support_v<1>(&prg_info.prg_markers_mask);
    prg_info.prg_markers_select = sdsl::select_support_mcl<1>(&prg_info.prg_markers_mask);

    prg_info.markers_mask_count_set_bits =
            prg_info.prg_markers_rank(prg_info.prg_markers_mask.size());

    prg_info.bwt_markers_mask = generate_bwt_markers_mask(prg_info.fm_index);

    prg_info.dna_bwt_masks = generate_bwt_masks(prg_info.fm_index, parameters);
    prg_info.rank_bwt_a = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_a);
    prg_info.rank_bwt_c = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_c);
    prg_info.rank_bwt_g = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_g);
    prg_info.rank_bwt_t = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_t);

    prg_info.num_variant_sites = prg_info.coverage_graph.bubble_map.size();
    return prg_info;
}
