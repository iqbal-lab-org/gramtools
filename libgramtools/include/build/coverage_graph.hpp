#ifndef COV_GRAPH_HPP
#define COV_GRAPH_HPP

#include "load_PRG_string.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <stack>

using namespace gram;
using seqPos = int32_t;

/**
 * The building blocks of a coverage graph
 * Contain sequence, site & allele ID, coverage array
 */
class coverage_Node {
    using covG_ptr = boost::shared_ptr<coverage_Node>;
public:
    coverage_Node() : sequence(""), site_ID(0), allele_ID(0), coverage(), pos(0), is_site_boundary{false} { ; };

    coverage_Node(seqPos pos) : sequence(""), site_ID(0), allele_ID(0), coverage(), pos(pos), is_site_boundary{false} { ; };

    coverage_Node(std::string const seq, int const pos, int const site_ID = 0, int const allele_ID = 0) :
            sequence(seq), pos(pos), site_ID(site_ID), allele_ID(allele_ID), is_site_boundary(false) { ; };

    bool is_boundary() { return is_site_boundary; };

    /**
     * Compare pointers to `coverage_Node`; used in topological ordering (lastmost sequence position first)
     * Equivalence in a set is defined using this, so we also test whether the pointers are the same objects.
     */
    friend bool operator>(const covG_ptr &lhs, const covG_ptr &rhs);

    bool has_sequence() { return sequence.size() != 0; };

    /*
     * Getters
     */
    int get_pos() { return pos; };

    /*
     * Setters
     */
    void set_pos(int pos) { this->pos = pos; };

    void add_sequence(std::string const &new_seq) { sequence += new_seq; };

    void add_edge(covG_ptr const target) { next.emplace_back(target); };

    void mark_as_boundary() { is_site_boundary = true; };

private:
    std::string sequence;
    int site_ID;
    int allele_ID;
    seqPos pos;
    std::vector<uint64_t> coverage;
    bool is_site_boundary;
    std::vector<covG_ptr> prev;
    std::vector<covG_ptr> next;
};

using covG_ptr = boost::shared_ptr<coverage_Node>;
using marker_to_node = std::unordered_map<Marker, covG_ptr>;

/**
* This class implements a DAG of `coverage_Node`s.
 * It is used to record coverage & to perform genotyping in gramtools.
**/
class coverage_Graph {
public:

    covG_ptr root;

    /**
     * Build a coverage graph from a PRG String int vector.
     */
    coverage_Graph(PRG_String const &vec_in);

    /** Maps the start of a local bubble, to its end.
     * Children nodes appear before parent nodes.
     */
    std::map<covG_ptr, covG_ptr, std::greater<covG_ptr> > bubble_map;

    parental_map par_map;
};

enum class marker_type {
    sequence, site_entry, allele_end, site_end
};

/**
 * A class in charge of mechanics of building coverage graph
 * It is designed for use by DEVELOPER only
 *
 */
class cov_Graph_Builder {
public:
    /*
     * FOR USER: All that's needed to make and run the object is under here
     */
    cov_Graph_Builder(PRG_String const &prg_string);

    void run();

    /*
     * 'Internal' functions
     */
    cov_Graph_Builder() {
        //std::cout << "WARNING: The builder needs a PRG_String object as constructor parameter to work properly.";
        ;
    };

    void make_root(); //* Start state: set up the globals such as `cur_Node` & `backWire`
    void make_sink(); //* End state: final wiring & pointers to null
    /**
     * Function call dispatcher based on marker at position @param pos.
     * Called once per element in `PRG_String`.
     * @param pos index into the `PRG_String`
     */
    void process_marker(uint32_t const &pos);
    void add_sequence(Marker const &m);
    marker_type find_marker_type(uint32_t const& pos);
    void enter_site(Marker const &m);
    void end_allele(Marker const &m);
    void exit_site(Marker const &m);
    /**
     * A convenience function for reaching the end of an allele: called by both
     * `end_allele` & `exit_site`
     */
    covG_ptr reach_allele_end(Marker const &m);
    void wire(covG_ptr const &target); // Build edges, 1 or 2, depending on whether `cur_Node` contains sequence.

    /*
     * variables & data structures
     */
    marker_vec linear_prg;
    std::unordered_map<Marker, int> end_positions;

    covG_ptr backWire; // Pointer to the most recent node needing edge building
    covG_ptr cur_Node;
    seqPos cur_pos; // For assigning position to nodes
    VariantLocus cur_Locus; // For building parental map

    marker_to_node bubble_starts;
    marker_to_node bubble_ends;

    // These will get transferred to the coverage_Graph
    covG_ptr root;
    std::map<covG_ptr, covG_ptr, std::greater<covG_ptr> > bubble_map;
    parental_map par_map;
};

#endif //COV_GRAPH_HPP
