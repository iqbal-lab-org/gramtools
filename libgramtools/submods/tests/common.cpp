#include "tests/common.hpp"
#include "genotype/quasimap/quasimap.hpp"

using namespace gram;

SitePbCoverage collect_coverage(coverage_Graph const& cov_graph, prg_positions positions){
    SitePbCoverage result(positions.size());
    covG_ptr accessed_node;
    std::size_t index{0};

    for (auto& pos : positions) {
        accessed_node = cov_graph.random_access[pos].node;
        result[index] = accessed_node->get_coverage();
        index++;
    }
    return result;
}

covG_ptrPair get_bubble_nodes(covG_ptr_map bubble_map, Marker site_ID){
    ensure_is_site_marker(site_ID);
    for (auto const& bubble : bubble_map){
        if (bubble.first->get_site_ID() == site_ID) return bubble;
    }
    throw std::invalid_argument("The provided site ID was not found in the map of PRG bubbles.");
}

void prg_setup::internal_setup(marker_vec encoded_prg, Sequences kmers){
    size_t kmer_size = kmers.front().size();
    for (auto const& kmer : kmers) assert(kmer_size == kmer.size());

    // TODO: the calls to rank_support setup in `generate_prg_info` get somehow lost when leaving its scope
    // and we need to call `init_support`, or rank_support again, in this scope for it to work
    prg_info = generate_prg_info(encoded_prg);

    sdsl::util::init_support(prg_info.rank_bwt_a, &prg_info.dna_bwt_masks.mask_a);
    sdsl::util::init_support(prg_info.rank_bwt_c, &prg_info.dna_bwt_masks.mask_c);
    sdsl::util::init_support(prg_info.rank_bwt_g, &prg_info.dna_bwt_masks.mask_g);
    sdsl::util::init_support(prg_info.rank_bwt_t, &prg_info.dna_bwt_masks.mask_t);

    sdsl::util::init_support(prg_info.prg_markers_rank, &prg_info.prg_markers_mask);
    sdsl::util::init_support(prg_info.prg_markers_select, &prg_info.prg_markers_mask);

    coverage = coverage::generate::empty_structure(prg_info);

    parameters.kmers_size = kmer_size;
    kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);
}


void prg_setup::quasimap_reads( GenomicRead_vector const& reads ){
    read_stats.compute_base_error_rate(reads);
    for (auto const& read : reads) {
        auto sequence = encode_dna_bases(read.seq);
        gram::quasimap_read(sequence, coverage, kmer_index, prg_info, parameters );
    }
    read_stats.compute_coverage_depth(coverage, prg_info.coverage_graph.par_map);
}
