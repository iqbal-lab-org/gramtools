/**
 * Produce graphviz dot file suitable for visualising the prg structure
 * and sequences.
 * Site entries/exits are labeled with the index they appear in in the PRG (and
 * the jvcf), and edges are labeled with the haplogroup of the allele series.
 */
#include <iostream>
#include <map>
#include <optional>
#include <regex>
#include <string>
#include <vector>

#include "genotype/infer/output_specs/segment_tracker.hpp"
#include "prg/coverage_graph.hpp"
#include "submod_resources.hpp"

using gram::genotype::SegmentTracker;
using gram::submods::covG_ptrPair;

void usage(const char* argv[]) {
  std::cout << "Usage: " << argv[0] << " coverage_graph region [coords_file]"
            << std::endl;
  std::cout << "\t coverage_graph: produced by gramtools build." << std::endl;
  std::cout << "\t region: description of subgraph to extract." << std::endl;
  std::cout << "\t\t region must be of form: 'chrom:start-stop' (for a genomic "
               "region) \n \t\t or "
               "'start-stop' (for site indices in the prg; 0-based, inclusive)."
            << std::endl;
  std::cout << "\t coords_file: produced by gramtools build. required if "
               "region is a genomic region."
            << std::endl;
  exit(1);
}

struct MatchRegion {
  std::string chrom;
  std::size_t start_location, end_location;
};

std::optional<MatchRegion> regexp_match_region(
    std::string const& region_string) {
  std::regex const genomic_region_regex("(.+):([0-9]+)-([0-9]+)");
  std::regex const site_index_region_regex("([0-9]+)-([0-9]+)");
  std::smatch pieces_match;
  std::size_t start_location, end_location;

  if (std::regex_match(region_string, pieces_match, genomic_region_regex)) {
    std::string chrom = pieces_match[1];
    start_location = stoi(pieces_match[2]);
    end_location = stoi(pieces_match[3]);
    return MatchRegion{chrom, start_location, end_location};
  } else if (std::regex_match(region_string, pieces_match,
                              site_index_region_regex)) {
    start_location = stoi(pieces_match[1]);
    end_location = stoi(pieces_match[2]);
    return MatchRegion{"", start_location, end_location};
  } else
    return {};
}

bool is_in_node(std::size_t query_pos, covG_ptr const& node,
                SegmentTracker& tracker) {
  auto const node_start = tracker.get_relative_pos(node->get_pos());
  auto const seq_size = node->get_sequence_size();
  auto node_stop = node_start;
  if (seq_size > 0) node_stop += seq_size - 1;
  return (query_pos >= node_start && query_pos <= node_stop);
}

std::optional<covG_ptrPair> find_nodes_by_genomic_region(
    MatchRegion const& match_region, covG_ptr const& root_node,
    SegmentTracker& tracker) {
  covG_ptr start_node, stop_node;
  std::size_t focal_position = match_region.start_location;
  covG_ptr cur_node = root_node;
  std::string chrom;
  while (cur_node->get_num_edges() > 0) {
    cur_node = cur_node->get_edges()[0];
    // TODO: below check is necessary because the sink node
    // is currently set to have the largest prg position + 1, triggering tracker
    // query to fail when no match found in the graph
    if (cur_node->get_num_edges() == 0) break;
    chrom = tracker.get_ID(cur_node->get_pos());
    if (chrom != match_region.chrom) continue;
    if (is_in_node(focal_position, cur_node, tracker)) {
      if (focal_position == match_region.start_location) {
        focal_position = match_region.end_location;
        start_node = cur_node;
      } else
        return covG_ptrPair{start_node, cur_node};
    }
  }
  return {};
}

covG_ptrPair find_nodes_by_site_idx(std::size_t const start_index,
                                    std::size_t const end_index,
                                    coverage_Graph const& graph) {
  auto num_var_sites = graph.bubble_map.size();
  if (start_index >= num_var_sites || end_index >= num_var_sites) {
    std::cout << "Error: there are only " + std::to_string(num_var_sites) +
                     " variant sites in the prg.\n";
    exit(1);
  }
  covG_ptr start_node, stop_node;
  auto start_pair = gram::submods::get_bubble_nodes(
      graph.bubble_map, index_to_siteID(start_index));
  start_node = start_pair.first;
  if (start_index != end_index)
    stop_node = gram::submods::get_bubble_nodes(graph.bubble_map,
                                                index_to_siteID(end_index))
                    .second;
  else
    stop_node = start_pair.second;
  return {start_node, stop_node};
}

int main(int argc, const char* argv[]) {
  if (argc < 3 || argc > 4) usage(argv);

  // Argument parsing and validation
  auto matched_region = regexp_match_region(argv[2]);
  if (!matched_region) {
    std::cout << "Error: invalid search region " << argv[2] << std::endl;
    usage(argv);
  }

  auto const start_location = matched_region->start_location;
  auto const end_location = matched_region->end_location;
  if (start_location > end_location) {
    std::cout << "Error: start must be <= stop \n";
    usage(argv);
  }

  SegmentTracker tracker;
  if (matched_region->chrom.size() != 0) {
    if (argc != 4) {
      std::cout << "Error: missing coords_file" << std::endl;
      usage(argv);
    }
    std::ifstream ifs{argv[3]};
    if (!ifs.good()) {
      std::cout << "Error: could not open " << argv[3] << std::endl;
      usage(argv);
    }
    tracker = SegmentTracker(ifs);
  }

  std::ifstream ifs{argv[1]};
  if (!ifs.good()) {
    std::cout << "Error: could not open " << argv[1] << std::endl;
    usage(argv);
  }
  boost::archive::binary_iarchive ia{ifs};
  coverage_Graph graph;
  ia >> graph;

  covG_ptrPair node_pair;
  if (matched_region->chrom.size() != 0) {
    auto result = find_nodes_by_genomic_region(matched_region.value(),
                                               graph.root, tracker);
    if (!result) {
      std::cout << "Error: could not find nodes spanning " << argv[2] << "\n";
      exit(1);
    }
    node_pair = result.value();
  } else {
    node_pair = find_nodes_by_site_idx(start_location, end_location, graph);
  }

  // Write subgraph by visiting all nodes inside the node pair
  auto start_node = node_pair.first;
  auto stop_node = node_pair.second;
  std::size_t new_idx{0};
  std::map<covG_ptr, std::size_t> node_ids{{start_node, new_idx++}};
  std::vector<covG_ptr> to_visit{start_node};

  covG_ptr cur_node;
  std::string nodes, edges, source, target;

  while (!to_visit.empty()) {
    cur_node = to_visit.back();
    to_visit.pop_back();
    if (cur_node->get_num_edges() == 0) continue;
    // Write node
    source = std::to_string(node_ids.at(cur_node));
    nodes.append(source);
    nodes.append(" [label=");
    if (cur_node->is_bubble_start() || cur_node->is_bubble_end())
      nodes.append(std::to_string(siteID_to_index(cur_node->get_site_ID())));
    else {
      auto node_seq = cur_node->get_sequence();
      if (cur_node->is_in_bubble()) {
        if (node_seq.empty()) node_seq = "\"\"";
        nodes.append(node_seq);
      } else {
        nodes.append("\"invariant ");
        nodes.append(std::to_string(node_seq.size()));
        nodes.append("bp\"");
      }
    }
    nodes.append("];\n");

    if (cur_node == stop_node) continue;

    std::size_t hapg{0};
    for (auto const& next_node : cur_node->get_edges()) {
      if (node_ids.find(next_node) == node_ids.end()) {
        node_ids.insert({next_node, new_idx++});
        to_visit.push_back(next_node);
      }

      // Write edge
      target = std::to_string(node_ids.at(next_node));
      edges.append(source);
      edges.append("->");
      edges.append(target);
      if (cur_node->is_bubble_start()) {
        edges.append(" [label=");
        edges.append(std::to_string(hapg));
        edges.append("]");
      }
      edges.append(";\n");
      hapg++;
    }
  }

  std::cout << "digraph \"gramtools_subgraph\" {\n" << nodes << edges << "\n}";
}
