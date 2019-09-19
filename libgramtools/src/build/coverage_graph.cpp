#include "build/coverage_graph.hpp"

coverage_Graph::coverage_Graph(PRG_String const &vec_in) {
    auto built_graph = cov_Graph_Builder(vec_in);
    root = built_graph.root;
    bubble_map = std::move(built_graph.bubble_map);
    par_map = std::move(built_graph.par_map);
    random_access = std::move(built_graph.random_access);
}

cov_Graph_Builder::cov_Graph_Builder(PRG_String const& prg_string) {
    linear_prg = prg_string.get_PRG_string();
    random_access = access_vec(linear_prg.size(), node_access());
    end_positions = prg_string.get_end_positions();
    make_root();
    cur_Locus = std::make_pair(0, 0); // Meaning: there is currently no current Locus.
}

void cov_Graph_Builder::run() {
    for(uint32_t i = 0; i < linear_prg.size(); ++i) process_marker(i);
    make_sink();
}

void cov_Graph_Builder::make_root() {
    cur_pos = -1;
   root = boost::make_shared<coverage_Node>(coverage_Node(cur_pos));
   backWire = root;
   cur_pos++;
   cur_Node = boost::make_shared<coverage_Node>(coverage_Node(cur_pos));
}

void cov_Graph_Builder::make_sink() {
    auto sink = boost::make_shared<coverage_Node>(coverage_Node(cur_pos + 1));
    wire(sink);
    cur_Node = nullptr;
    backWire = nullptr;
}

void cov_Graph_Builder::process_marker(uint32_t const &pos) {
    Marker m = linear_prg[pos];
    marker_type t = find_marker_type(pos);

    switch(t){
        case marker_type::sequence:
            add_sequence(m);
            break;
        case marker_type::site_entry:
            enter_site(m);
            break;
        case marker_type::allele_end:
            end_allele(m);
            break;
        case marker_type::site_end:
            exit_site(m);
    }

    // Set up random access
    auto seq_size = cur_Node->get_sequence_size();
    if (seq_size <= 1) // Will include all site entry and exit nodes, and sequence nodes with a single character
        random_access[pos] = node_access{cur_Node, 0};
    else random_access[pos] = node_access{cur_Node, seq_size - 1};
}

marker_type cov_Graph_Builder::find_marker_type(uint32_t const& pos) {
    auto const& m = linear_prg[pos];
    if (m <= 4) return marker_type::sequence; // Note: the `PRG_String` constructor code must make sure that m is > 0.

    // After passing through `PRG_String` constructor code, the only time a marker is odd is when it signals a site entry.
    if (m % 2 == 1) return marker_type::site_entry;

    // Find if the allele marker, which must be in the map, is at the end position or not
    auto& end_pos = end_positions.at(m);
    assert(pos <= end_pos); // Sanity check: the allele marker must occur before, or at, the end position of the site
    if (pos < end_pos) return marker_type::allele_end;

    return marker_type::site_end;
}

void cov_Graph_Builder::add_sequence(Marker const &m) {
    std::string c = decode_dna_base(m); // Note: implicit conversion of m from uint32_t to uint8_t; this is OK because 0 < m < 5.
    cur_Node->add_sequence(c);
    cur_pos++;
}

void cov_Graph_Builder::enter_site(Marker const &m) {

    auto site_entry = boost::make_shared<coverage_Node>(coverage_Node("", cur_pos, m, 0));
    site_entry->mark_as_boundary();
    wire(site_entry);

    // Update the global pointers
    cur_Node = boost::make_shared<coverage_Node>(coverage_Node("", cur_pos, m, 1));
    backWire = site_entry;

    // Make & register a new bubble
    auto site_exit = boost::make_shared<coverage_Node>(coverage_Node("", cur_pos, m, 0));
    site_exit->mark_as_boundary();
    bubble_map.insert(std::make_pair(site_entry,site_exit));
    bubble_starts.insert(std::make_pair(m, site_entry));
    bubble_ends.insert(std::make_pair(m, site_exit));

    // Update the parent map & the current Locus
    if (cur_Locus.first != 0) {
        assert(par_map.find(m) == par_map.end()); // The marker should not be in there already
        par_map.insert(std::make_pair(m, cur_Locus));
    }
    cur_Locus = std::make_pair(m, 1);
}

void cov_Graph_Builder::end_allele(Marker const &m) {
    auto site_ID = m - 1;
    auto site_exit = reach_allele_end(m);

    auto& allele_ID = cur_Locus.second;
    // Reset node and position to the site start node
    auto site_entry = bubble_starts.at(site_ID);
    backWire = site_entry;
    cur_pos = site_entry->get_pos();

    // Update to the next allele
    allele_ID++; // increment allele ID.
    cur_Node = boost::make_shared<coverage_Node>(coverage_Node("", cur_pos, site_ID, allele_ID));
}

void cov_Graph_Builder::exit_site(Marker const &m) {
    auto site_ID = m - 1;
    auto site_exit = reach_allele_end(m);

    int allele_ID;
    // Update the current Locus
    try {
        cur_Locus = par_map.at(site_ID);
        allele_ID = cur_Locus.second;
    }
    catch(std::out_of_range&e){ // Means we were in a level 1 site; we will no longer be in a site
        cur_Locus = std::make_pair(0,0);
        allele_ID = 0;
    }

    backWire = site_exit;
    cur_pos = site_exit->get_pos(); // Take the largest allele pos as the new current pos.
    cur_Node = boost::make_shared<coverage_Node>(coverage_Node("", cur_pos, site_ID, allele_ID));
}

covG_ptr cov_Graph_Builder::reach_allele_end(Marker const& m){
    // Make sure we are tracking the right site
    auto site_ID = m - 1;
    assert(cur_Locus.first == site_ID);

    auto site_exit = bubble_ends.at(site_ID);
    wire(site_exit);

    // Update the exit's pos if it is smaller than this allele's
    if (site_exit->get_pos() < cur_pos) site_exit->set_pos(cur_pos);

    return site_exit;
}


void cov_Graph_Builder::wire(covG_ptr const &target) {
    if (cur_Node->has_sequence()){
        backWire->add_edge(cur_Node);
        cur_Node->add_edge(target);
    }
    else backWire->add_edge(target);
}

bool operator>(const covG_ptr &lhs, const covG_ptr &rhs) {
    if (lhs->pos == rhs->pos) {
        return lhs.get() > rhs.get();
    } else return lhs->pos > rhs->pos;
}
