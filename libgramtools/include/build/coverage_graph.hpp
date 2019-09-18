#include "load_PRG_string.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <stack>

using namespace gram;

enum class marker_type{sequence, site_entry, allele_end, site_end};

/**
 * The building blocks of a coverage graph
 * Contain sequence, site & allele ID, coverage array
 */
class coverage_Node{
    using covG_ptr = boost::shared_ptr<coverage_Node>;
public:
    coverage_Node() : sequence(""), site_ID(0), allele_ID(0), coverage(), is_site_boundary{false}{;};
    coverage_Node(std::string const seq, int const pos, int const site_ID = 0, int const allele_ID = 0);
    bool is_boundary() {return is_site_boundary;};

    /**
     * Compare pointers to `coverage_Node`; used in topological ordering (lastmost sequence position first)
     * Equivalence in a set is defined using this, so we also test whether the pointers are the same objects.
     */
    friend bool operator > (const covG_ptr& lhs, const covG_ptr& rhs);

    bool has_sequence(){return sequence.size() != 0;};
    /*
     * Getters
     */
    int get_pos(){return pos;};

    /*
     * Setters
     */
    void set_pos(int pos){this->pos = pos;};
    void add_sequence(std::string const& new_seq){sequence += new_seq;};
    void add_edge(covG_ptr const target) {next.emplace_back(target);};
    void mark_as_boundary(){is_site_boundary = true ;};

private:
    std::string sequence;
    int site_ID;
    int allele_ID;
    int pos;
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
class coverage_Graph{
public:

    covG_ptr root;

    /**
     * Build a coverage graph from a PRG String int vector.
     */
    coverage_Graph(PRG_String const& vec_in);

    /** Maps the start of a local bubble, to its end.
     * Children nodes appear before parent nodes.
     */
    std::map<covG_ptr,covG_ptr, std::greater<covG_ptr> > bubble_map;

    parental_map par_map;
};

/**
 * A class in charge of mechanics of building coverage graph
 * It is designed for DEVELOPER only
 *
 */
class cov_Graph_Builder{
    cov_Graph_Builder(PRG_String const& prg_string);

    /*
     * functions
     */
    void make_root();
    void make_sink();
    /**
     * Function call dispatcher based on marker at position @param pos.
     * Called once per element in `PRG_String`.
     * @param pos index into the `PRG_String`
     */
    void process_marker(uint32_t const& pos);
    void add_sequence(Marker const& m);
    marker_type find_marker_type(Marker const& m);
    void enter_site(Marker const& m);
    void end_allele(Marker const& m);
    void exit_site(Marker const& m);
    void wire(covG_ptr const& target);

    /*
     * variables & data structures
     */
    PRG_String const* prg_string; // Pointer to the constructor parameter
    covG_ptr backWire;
    covG_ptr cur_Node;
    uint32_t cur_pos; // For assigning position to nodes
    VariantLocus cur_Locus; // For building parental map

    marker_to_node bubble_starts;
    marker_to_node bubble_ends;

    // These will get transferred to the coverage_Graph
    covG_ptr root;
    std::map<covG_ptr,covG_ptr, std::greater<covG_ptr> > bubble_map;
    parental_map par_map;
};
