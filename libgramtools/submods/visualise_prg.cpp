/**
 * Produce graphviz dot file suitable for visualising the prg structure
 * and sequences.
 */
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "prg/coverage_graph.hpp"

#include "submod_resources.hpp"

void usage(const char* argv[]) {
  std::cout << "Usage: " << argv[0] << " coverage_graph start stop outfile"
            << std::endl;
  std::cout << "\t coverage_graph: produced by gramtools build." << std::endl;
  std::cout << "\t start: index of first site to visualise" << std::endl;
  std::cout << "\t stop: index of last site to visualise. can be same as start."
            << std::endl;
  exit(1);
}

int main(int argc, const char* argv[]) {
  if (argc != 5) usage(argv);

  // Argument parsing and validation
  int start_idx = atoi(argv[2]);
  int stop_idx = atoi(argv[3]);
  if (start_idx < 0 || stop_idx < 0 || (stop_idx < start_idx)) {
    std::cout << "Error: start and stop must be >0 and stop >= start \n";
    usage(argv);
  }

  std::string ofprefix{argv[4]};
  std::ofstream ofs{ofprefix + ".gv"};
  if (!ofs.good()) {
    std::cout << "Error: could not open " << argv[4] << std::endl;
    usage(argv);
  }

  std::ifstream ifs{argv[1]};
  if (!ifs.good()) {
    std::cout << "Error: could not open " << argv[1] << std::endl;
    usage(argv);
  }
  boost::archive::binary_iarchive ia{ifs};
  coverage_Graph graph;
  ia >> graph;

  auto num_var_sites = graph.bubble_map.size();
  if (start_idx >= num_var_sites || stop_idx >= num_var_sites) {
    std::cout << "Error: there are only " + std::to_string(num_var_sites) +
                     " variant sites in the prg.\n";
    usage(argv);
  }

  covG_ptr start_node, stop_node;
  auto start_pair = gram::submods::get_bubble_nodes(graph.bubble_map,
                                                    index_to_siteID(start_idx));
  start_node = start_pair.first;
  if (start_idx != stop_idx)
    stop_node = gram::submods::get_bubble_nodes(graph.bubble_map,
                                                index_to_siteID(stop_idx))
                    .second;
  else
    stop_node = start_pair.second;

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
      nodes.append("\"\"");
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
      edges.append(";\n");
    }
  }

  ofs << "digraph \"" << ofprefix << "\" {\n" << nodes << edges << "\n}";
}
