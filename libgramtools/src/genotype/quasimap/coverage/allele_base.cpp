#include <cassert>
#include <fstream>
#include <vector>

#include "genotype/quasimap/coverage/allele_base.hpp"


using namespace gram;
using namespace gram::coverage::per_base;


SitesAlleleBaseCoverage gram::coverage::generate::allele_base_non_nested(const PRG_Info &prg_info) {
    // If graph is nested, this data structure cannot be populated correctly, so return it empty by convention
    if (prg_info.coverage_graph.is_nested) return SitesAlleleBaseCoverage{};

    uint64_t number_of_variant_sites = prg_info.num_variant_sites;
    SitesAlleleBaseCoverage allele_base_coverage(number_of_variant_sites);

    Marker site_ID;
    const auto min_boundary_marker = 5;

    for (auto const& bubble_entry : prg_info.coverage_graph.bubble_map){
        site_ID = bubble_entry.first->get_site_ID();
        auto site_ID_corresp_index = (site_ID - min_boundary_marker) / 2;
        SitePbCoverage& referent = allele_base_coverage.at(site_ID_corresp_index);

        for (auto const& allele_node : bubble_entry.first->get_edges()){
            assert(allele_node->is_in_bubble());
            // Add as many 0's in the allele as there are bases in the node
           referent.emplace_back(PerBaseCoverage(allele_node->get_coverage()));
        }
    }
    return allele_base_coverage;
}

void coverage::record::allele_base(PRG_Info const &prg_info, const SearchStates &search_states,
                                   const uint64_t &read_length) {
    PbCovRecorder record_it{prg_info, search_states, read_length};
}

/**
 * String serialise the base coverages for one allele.
 */
std::string dump_allele(const PerBaseCoverage &allele) {
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
std::string dump_site(const SitePbCoverage &site) {
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
        cur_Node(start_point.node), traversed_loci(traversed_loci), bases_remaining(read_size), first_node(true),
        end_pos(0){

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
    end_pos = 0;
    std::size_t seq_size = cur_Node->get_sequence_size();
    if (seq_size > 0) end_pos = std::min(seq_size - 1, start_pos + bases_remaining - 1);
}

void Traverser::choose_allele() {
    auto traversed_locus = traversed_loci[traversed_index];
    auto site_id{traversed_locus.first}, allele_id{traversed_locus.second};
    auto next_node = cur_Node->get_edges()[allele_id - 1];

    // Check site & allele consistency
    if (next_node->has_sequence()){
        assert(next_node->get_site_ID() == site_id && next_node->get_allele_ID() == allele_id);
    }

    cur_Node = next_node;
}


PbCovRecorder::PbCovRecorder(const PRG_Info &prg_info, SearchStates const &search_states,
                             std::size_t read_size)
: prg_info(&prg_info), read_size(read_size){
        for (auto const& search_state : search_states) process_SearchState(search_state);
        write_coverage_from_dummy_nodes();
}

void PbCovRecorder::write_coverage_from_dummy_nodes(){
    covG_ptr cov_node;
    node_coordinates to_increment;
   for (auto const& element : cov_mapping){ // Go through each dummy node
       cov_node = element.first;
       to_increment = element.second.get_coordinates();
       PerBaseCoverage& cur_coverage = cov_node->get_ref_to_coverage(); // Modifiable in place
       for (auto i = to_increment.first; i <= to_increment.second; i++) {
           if (cur_coverage[i] == UINT16_MAX) continue;
#pragma omp atomic
           cur_coverage[i]++;
       }
   }
}


void PbCovRecorder::process_SearchState(SearchState const& ss){
    bool first{true};
    Traverser t;

    for (auto occurrence = ss.sa_interval.first; occurrence <= ss.sa_interval.second; occurrence++){
        auto coordinate = prg_info->fm_index[occurrence];
        auto access_point = prg_info->coverage_graph.random_access[coordinate];
        t = {access_point, ss.traversed_path, read_size};

        // Record a full traversal starting at the first mapping instance
        if (first){
            first = false;
            record_full_traversal(t);
        }
        // Record coverage at the alternative start points in the SearchState, if any.
        // This occurs if the suffix interval has size more than one.
        // In which case we only want need to process the first node in the traversal
        else {
            auto cur_Node = t.next_Node().value();
            auto coordinates = t.get_node_coordinates();
            process_Node(cur_Node, coordinates.first, coordinates.second);
        }
    }
}

void PbCovRecorder::record_full_traversal(Traverser& t){
    auto cur_Node = t.next_Node();
    auto coordinates = t.get_node_coordinates();
    while(bool(cur_Node)){
        process_Node(cur_Node.value(), coordinates.first, coordinates.second);
        cur_Node = t.next_Node();
        coordinates = t.get_node_coordinates();
    }
}

void PbCovRecorder::process_Node(covG_ptr cov_node, node_coordinate start_pos, node_coordinate end_pos){
    if (! cov_node->has_sequence()) return; // Skips double site entries, where `cov_node` is a no-sequence bubble entry
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
