#include <cassert>
#include <fstream>
#include <vector>

#include "quasimap/coverage/allele_base.hpp"


using namespace gram;
using namespace gram::coverage::per_base;


SitesAlleleBaseCoverage gram::coverage::generate::allele_base_structure(const PRG_Info &prg_info) {
    uint64_t number_of_variant_sites = prg_info.num_variant_sites;
    SitesAlleleBaseCoverage allele_base_coverage(number_of_variant_sites);

    const auto min_boundary_marker = 5;

    uint64_t allele_size = 0;
    Marker last_marker = 0;

    // Traverse the sites mask, in order to identify alleles.
    for (const auto &mask_value: prg_info.sites_mask) {
        auto within_allele = mask_value != 0;
        if (within_allele) {
            allele_size += 1;
            last_marker = mask_value;
            continue;
        }

        auto no_allele_to_flush = allele_size == 0;
        if (no_allele_to_flush)
            continue;

        // Store room aside for the allele
        BaseCoverage bases(allele_size);
        uint64_t variant_site_cover_index = (last_marker - min_boundary_marker) / 2;
        allele_base_coverage.at(variant_site_cover_index).emplace_back(bases);
        allele_size = 0;
    }
    return allele_base_coverage;
}


uint64_t gram::allele_start_offset_index(const uint64_t within_allele_prg_index, const PRG_Info &prg_info) {
    uint64_t number_markers_before = prg_info.prg_markers_rank(within_allele_prg_index);
    // Rank operation gets index of nearest left marker in prg, marking allele's start.
    uint64_t marker_index = prg_info.prg_markers_select(number_markers_before);
    uint64_t offset = within_allele_prg_index - marker_index - 1;

    return offset;
}


uint64_t gram::set_site_base_coverage(Coverage &coverage,
                                      SitesCoverageBoundaries &sites_coverage_boundaries,
                                      const VariantLocus &path_element,
                                      const uint64_t allele_coverage_offset,
                                      const uint64_t max_bases_to_set) {
    // Extract the variant site of interest using the variant site marker number.
    auto marker = path_element.first;
    auto min_boundary_marker = 5;
    auto variant_site_coverage_index = (marker - min_boundary_marker) / 2;
    auto &site_coverage = coverage.allele_base_coverage.at(variant_site_coverage_index);

    // Extract the allele of interest using the allele id.
    auto allele_id = path_element.second;
    auto allele_coverage_index = allele_id - 1;
    auto &allele_coverage = site_coverage.at(allele_coverage_index);

    // Now: which bases inside the allele are covered by the read?
    // If `index_end_boundary` gets set to `allele_coverage_offset+max_bases_to_set`, the read ends before the allele's end.
    uint64_t index_end_boundary = std::min(allele_coverage_offset + max_bases_to_set, allele_coverage.size());
    assert(index_end_boundary >= allele_coverage_offset);
    uint64_t count_bases_consumed = index_end_boundary - allele_coverage_offset;

    uint64_t index_start_boundary = allele_coverage_offset;
    bool site_seen_previously = sites_coverage_boundaries.find(path_element) != sites_coverage_boundaries.end();
    // If we have already mapped to this `VariantLocus` before, we allow only to map from the end of the previous mapping onwards.
    // TODO: what if we map the end of the allele previously and now map the whole allele? Here, we could not map the missed bases at the beginning...
    if (site_seen_previously)
        index_start_boundary = std::max(allele_coverage_offset, sites_coverage_boundaries[path_element]);
    sites_coverage_boundaries[path_element] = index_end_boundary; // Update the end_index mapped.

    // Actually increment the base counts between specified ranges.
    for (uint64_t i = index_start_boundary; i < index_end_boundary; ++i) {
        if (allele_coverage[i] == UINT16_MAX)
            continue;

#pragma omp atomic
        ++allele_coverage[i];
    }
    return count_bases_consumed;
}


std::pair<uint64_t, uint64_t> gram::site_marker_prg_indexes(const uint64_t &site_marker, const PRG_Info &prg_info) {
    auto alphabet_rank = prg_info.fm_index.char2comp[site_marker];
    auto first_sa_index = prg_info.fm_index.C[alphabet_rank];

    auto first_prg_index = prg_info.fm_index[first_sa_index];
    // Need to be sure we are dealing with a site marker so that we know we can look for its even counterpart in map
    assert(is_site_marker(site_marker));
    auto second_prg_index = prg_info.last_allele_positions.at(site_marker + 1);

    return std::make_pair(first_prg_index, second_prg_index);
}


/**
 * A read may map to a set of variant sites in different ways (by traversing different alleles).
 * Here we record base level coverage for one such traversal (ie a `SearchState`)
 * The details of this function are here only to deal with reads that:
 * * Start inside an allele
 * * End inside an allele
 * If not, we increment all bases inside each traversed allele- but only once!
 * @param sites_coverage_boundaries: hashes a `VariantLocus` to the last base in the allele with recorded coverage
 * If this base is the last base of the allele, allows to not record coverage of an allele more than once.
 */
void sa_index_allele_base_coverage(Coverage &coverage,
                                   SitesCoverageBoundaries &sites_coverage_boundaries,
                                   const uint64_t &sa_index,
                                   const uint64_t &read_length,
                                   const SearchState &search_state,
                                   const PRG_Info &prg_info) {
    uint64_t read_bases_consumed = 0;
    uint64_t last_site_marker = 0;
    std::pair<uint64_t, uint64_t> last_site_prg_start_end = std::make_pair(0, 0);
    std::pair<uint64_t, uint64_t> site_prg_start_end = std::make_pair(0, 0);
    auto path_it = search_state.traversed_path.rbegin();

    auto read_start_index = prg_info.fm_index[sa_index]; // Where the mapping instance starts in the prg.
    auto start_site_marker = prg_info.sites_mask[read_start_index];
    bool read_starts_within_site = start_site_marker != 0; // Are we inside a variant site?
    if (read_starts_within_site) {
        const auto &path_element = *path_it;
        auto site_marker = path_element.first;
        last_site_prg_start_end = site_marker_prg_indexes(site_marker, prg_info);

        auto allele_coverage_offset = allele_start_offset_index(read_start_index, prg_info);
        auto max_bases_to_set = read_length - read_bases_consumed;
        read_bases_consumed += set_site_base_coverage(coverage,
                                                      sites_coverage_boundaries,
                                                      path_element,
                                                      allele_coverage_offset,
                                                      max_bases_to_set);
        ++path_it;
    } else {
        // Fast-forward to next variant site. Just need to consume bases going up to there.
        const auto &path_element = *path_it;
        auto site_marker = path_element.first;
        site_prg_start_end = site_marker_prg_indexes(site_marker, prg_info);
        read_bases_consumed += site_prg_start_end.first - read_start_index;
    }

    auto last_path_it = search_state.traversed_path.rend();
    while (read_bases_consumed < read_length and path_it != last_path_it) {
        const auto &path_element = *path_it;
        auto site_marker = path_element.first;

        // Case: consume invariant bases between two variant sites
        if (last_site_prg_start_end.first != 0) {
            site_prg_start_end = site_marker_prg_indexes(site_marker, prg_info);
            read_bases_consumed += site_prg_start_end.first - last_site_prg_start_end.second - 1;
            if (read_bases_consumed > read_length)
                throw std::out_of_range(
                        "ERROR:Consumed more bases than the read length."
                        "Check that the site boundaries are properly detected.");
        }
        last_site_prg_start_end = site_prg_start_end;

        uint64_t allele_coverage_offset = 0;
        auto max_bases_to_set = read_length - read_bases_consumed;
        read_bases_consumed += set_site_base_coverage(coverage,
                                                      sites_coverage_boundaries,
                                                      path_element,
                                                      allele_coverage_offset,
                                                      max_bases_to_set);
        ++path_it;
    }
}


void coverage::record::allele_base(Coverage &coverage,
                                   const SearchStates &search_states,
                                   const uint64_t &read_length,
                                   const PRG_Info &prg_info) {
    SitesCoverageBoundaries sites_coverage_boundaries;

    /*
     class NodeImprint{
        cov_Node* referent;
        uint32_t start;
         uint32_t end;
         bool full; // if end-start + 1 == referent.coverage.size()
     }

     class pbcovRecorder{
         void generate_imprints(SearchState const& ss);
         void record_coverage(corresp);

         cov_graph* cov_graph;
         std::map<cov_Node*, NodeImprint> corresp;
         SearchStates input_ss;
     };
     */

    for (const auto &search_state: search_states) {
        if (search_state.traversed_path.empty())
            continue;

        auto first_sa_index = search_state.sa_interval.first;
        auto last_sa_index = search_state.sa_interval.second;
        // Record base-level coverage for each mapped instance of the read.
        // TODO: can you have a sa_index bigger than one for SearchState with non empty variant_site_path?
        for (auto sa_index = first_sa_index; sa_index <= last_sa_index; ++sa_index)
            sa_index_allele_base_coverage(coverage,
                                          sites_coverage_boundaries,
                                          sa_index,
                                          read_length,
                                          search_state,
                                          prg_info);
    }
}

/**
 * String serialise the base coverages for one allele.
 */
std::string dump_allele(const BaseCoverage &allele) {
    std::stringstream stream;
    stream << "[";
    auto i = 0;
    for (const auto &base_coverage: allele) {
        stream << (int) base_coverage;
        if (i++ < allele.size() - 1)
            stream << ",";
    }
    stream << "]";
    return stream.str();
}

/**
 * String serialise the alleles of a site.
 * @see dump_allele()
 */
std::string dump_site(const AlleleCoverage &site) {
    std::stringstream stream;
    auto i = 0;
    for (const auto &allele: site) {
        stream << dump_allele(allele);
        if (i++ < site.size() - 1)
            stream << ",";
    }
    return stream.str();
}

/**
 * String serialise all base-level coverages for all sites of the prg.
 * @see dump_site()
 */
std::string dump_sites(const SitesAlleleBaseCoverage &sites) {
    std::stringstream stream;
    auto i = 0;
    for (const auto &site: sites) {
        stream << "[";
        stream << dump_site(site);
        stream << "]";
        if (i++ < sites.size() - 1)
            stream << ",";
    }
    return stream.str();
}

std::string gram::dump_allele_base_coverage(const SitesAlleleBaseCoverage &sites) {
    std::stringstream stream;
    stream << "{\"allele_base_counts\":[";
    stream << dump_sites(sites);
    stream << "]}";
    return stream.str();
}

void coverage::dump::allele_base(const Coverage &coverage,
                                 const Parameters &parameters) {
    std::string json_string = dump_allele_base_coverage(coverage.allele_base_coverage);
    std::ofstream file;
    file.open(parameters.allele_base_coverage_fpath);
    file << json_string << std::endl;
}

DummyCovNode::DummyCovNode(node_coordinate start_pos, node_coordinate end_pos, std::size_t node_size)
: start_pos(start_pos), end_pos(end_pos), node_size(node_size), full(false){
    if (start_pos > end_pos) throw InconsistentCovNodeCoordinates("start_pos must not be greater than end_pos");
    if (start_pos >= node_size || end_pos >= node_size) throw \
        InconsistentCovNodeCoordinates("node_size must be greater than start_pos and end_pos");

    if(end_pos - start_pos == node_size - 1) full = true;
}

void DummyCovNode::extend_coordinates(node_coordinates coords) {
    if (coords.second >= node_size) throw InconsistentCovNodeCoordinates("end coordinate must be less than node_size");
    if (full) return;
    if (coords.first < start_pos) start_pos = coords.first;
    if (coords.second > end_pos){
        end_pos = coords.second;
    }
    if (end_pos - start_pos == node_size - 1) full = true;
}


Traverser::Traverser(node_access start_point, VariantSitePath traversed_loci, std::size_t read_size) :
        cur_Node(start_point.node), traversed_loci(traversed_loci), bases_remaining(read_size), first_node(true){

    // Will give the next site to choose in the graph
    traversed_index = traversed_loci.size();
    start_pos = start_point.offset;
}


std::optional<covG_ptr> Traverser::next_Node() {
    if (first_node){
        process_first_node();
        first_node = false;
        return cur_Node;
    }
    else if (bases_remaining == 0) {
        return {};
    }
    else {
        go_to_next_site();
        if (cur_Node == nullptr) return {};
        return cur_Node;
    }
}

void Traverser::process_first_node(){
    update_coordinates();
    if (!cur_Node->is_in_bubble()) go_to_next_site();
}

void Traverser::go_to_next_site(){
    start_pos = 0;
    // Skip invariants
    while(cur_Node->get_edges().size() == 1){
        if (bases_remaining <= 0) {
            cur_Node = nullptr;
            return;
        }
        move_past_single_edge_node();
        update_coordinates();
        if (cur_Node->is_in_bubble()) return; // Deals with exiting nested sites: we need to avoid skipping those
    }
    // Pick the right allelic node
    --traversed_index;
    choose_allele();

    update_coordinates();
}

void Traverser::update_coordinates(){
    assign_end_position();
    if (cur_Node->has_sequence()) bases_remaining -= (end_pos - start_pos + 1);
}

void Traverser::move_past_single_edge_node(){
    assert(cur_Node->get_edges().size() == 1);
    cur_Node = cur_Node->get_edges()[0];
}

void Traverser::assign_end_position() {
    std::size_t seq_size = cur_Node->get_sequence_size();
    end_pos = std::min(seq_size - 1, start_pos + bases_remaining - 1);
}

void Traverser::choose_allele() {
    auto traversed_locus = traversed_loci[traversed_index];
    auto site_id{traversed_locus.first}, allele_id{traversed_locus.second};
    auto next_node = cur_Node->get_edges()[allele_id - 1];
    assert(next_node->get_site() == site_id && next_node->get_allele() == allele_id);

    cur_Node = next_node;
}


PbCovRecorder::PbCovRecorder(PRG_Info& prg_info, SearchStates const& search_states,
       std::size_t read_size)
: prg_info(&prg_info), read_size(read_size){
        for (auto const& search_state : search_states) process_SearchState(search_state);
}

void PbCovRecorder::process_SearchState(SearchState const& ss){
    std::size_t sa_interval_size = ss.sa_interval.second - ss.sa_interval.first + 1;
    auto first_coordinate = prg_info->fm_index[ss.sa_interval.first];
    auto first_access_point = prg_info->coverage_graph.random_access[first_coordinate];
    Traverser t{first_access_point, ss.traversed_path, read_size};

    auto cur_Node = t.next_Node();
    auto coordinates = t.get_node_coordinates();
    while(bool(cur_Node)){
        process_Node(cur_Node.value(), coordinates.first, coordinates.second);
        cur_Node = t.next_Node();
        coordinates = t.get_node_coordinates();
    }
}

void PbCovRecorder::process_Node(covG_ptr cov_node, node_coordinate start_pos, node_coordinate end_pos){
    if (cov_mapping.find(cov_node) == cov_mapping.end()){
        std::size_t cov_node_size = cov_node->get_sequence_size();
        DummyCovNode new_dummy_cov_node{start_pos, end_pos, cov_node_size};
        cov_mapping.insert({cov_node, new_dummy_cov_node});
    }
    else{
        DummyCovNode& existing_dummy_cov_node = cov_mapping.at(cov_node);
        existing_dummy_cov_node.extend_coordinates(node_coordinates{start_pos, end_pos});
    }
}
