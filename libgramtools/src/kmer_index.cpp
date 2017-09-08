#include <memory>

#include "bidir_search_bwd.hpp"
#include "kmer_index.hpp"


using Kmer = std::vector<uint8_t>;
using Kmers = std::vector<Kmer>;

using SA_Interval = std::pair<uint64_t, uint64_t>;
using SA_Intervals = std::list<SA_Interval>;


struct SA_IntervalsCacheNode {
    std::unique_ptr<SA_Intervals> sa_intervals;
    std::unique_ptr<Sites> sites;
    uint8_t _dbg_base = 0;

    std::weak_ptr<SA_IntervalsCacheNode> parent;
    std::array<std::shared_ptr<SA_IntervalsCacheNode>, 4> children;
    bool all_children_have_sa_intervals = false;
};


void print_sa_intervals(const SA_Intervals &sa_intervals) {
    for (const auto &sa_interval: sa_intervals)
        std::cout << "(" << sa_interval.first << ", " << sa_interval.second << ")" << "   ";
}


std::string kmer_cache_to_str(const std::shared_ptr<SA_IntervalsCacheNode> root);


std::shared_ptr<SA_IntervalsCacheNode> construct_child_node(const std::shared_ptr<SA_IntervalsCacheNode> parent_node,
                                                            const uint8_t kmer_base,
                                                            const uint64_t max_alphabet_num,
                                                            const std::vector<int> &allele_mask,
                                                            const DNA_Rank &rank_all,
                                                            const FM_Index &fm_index) {
    auto child_node = std::make_shared<SA_IntervalsCacheNode>();
    child_node->parent = parent_node;
    child_node->_dbg_base = kmer_base;

    SA_Intervals sa_intervals = *(parent_node->sa_intervals);
    Sites sites = *(parent_node->sites);

    reduce_sa_intervals(kmer_base,
                        sa_intervals,
                        sites,
                        false, true,
                        allele_mask,
                        max_alphabet_num,
                        false,
                        rank_all,
                        fm_index);

    std::cout << (int) child_node->_dbg_base << ": size " << sa_intervals.size() << std::endl;
    print_sa_intervals(sa_intervals);
    std::cout << std::endl;

    child_node->sa_intervals = std::make_unique<SA_Intervals>(std::move(sa_intervals));
    child_node->sites = std::make_unique<Sites>(std::move(sites));

    return child_node;
}


void cache_kmer(const Kmer &kmer,
                std::shared_ptr<SA_IntervalsCacheNode> kmer_cache_root,
                const uint64_t max_alphabet_num,
                const std::vector<int> &allele_mask,
                const DNA_Rank &rank_all,
                const FM_Index &fm_index) {

    auto node = kmer_cache_root;

    for (const auto &kmer_base: kmer) {
        const int child_node_idx = (int) kmer_base - 1;
        auto &next_node = node->children[child_node_idx];
        if (next_node) {
            node = next_node;
            continue;
        }

        next_node = construct_child_node(node,
                                         kmer_base,
                                         max_alphabet_num,
                                         allele_mask,
                                         rank_all,
                                         fm_index);
        node = next_node;
    }
}


void generate_kmer_index(const Kmers &kmers,
                         const uint64_t max_alphabet_num,
                         const std::vector<int> &allele_mask,
                         const DNA_Rank &rank_all,
                         const FM_Index &fm_index) {

    auto kmer_cache_root = std::make_shared<SA_IntervalsCacheNode>();

    SA_Intervals sa_intervals = {{0, fm_index.size()}};
    kmer_cache_root->sa_intervals = std::make_unique<SA_Intervals>(std::move(sa_intervals));

    Sites sites = {Site()};
    kmer_cache_root->sites = std::make_unique<Sites>(std::move(sites));

    for (const auto &kmer: kmers)
        cache_kmer(kmer, kmer_cache_root,
                   max_alphabet_num,
                   allele_mask,
                   rank_all,
                   fm_index);

    auto cache_str = kmer_cache_to_str(kmer_cache_root);
    std::cout << cache_str << std::endl;
}


bool node_is_leaf(const std::shared_ptr<SA_IntervalsCacheNode> node) {
    for (auto i = 0; i < 4; ++i) {
        auto child = node->children[i];
        if (child) {
            return false;
        }
    }
    return true;
}


std::string kmer_cache_to_str(const std::shared_ptr<SA_IntervalsCacheNode> root) {
    std::stringstream s;
    if (node_is_leaf(root)) {
        s << (int) root->_dbg_base;
        return s.str();
    }

    s << (int) root->_dbg_base << "[ ";
    for (auto i = 0; i < 4; ++i) {
        auto node = root->children[i];

        if (!node) {
            s << ".";
            if (i != 3)
                s << "  ";
            continue;
        }

        auto rep = kmer_cache_to_str(node);
        s << rep;
        if (i != 3)
            s << "  ";
    }
    s << " ]";

    return s.str();
}


/*
 * TODO:
 * Investigate use of csa_bitcompressed<int_alphabet<> > equivalent for 4 char DNA bases as
 * first template parameter for cst.


using CST_Type = sdsl::cst_sct3<sdsl::csa_bitcompressed<sdsl::int_alphabet<> > >;

void traverse(CST_Type &cst);

template<class t_cst, class t_pat=typename t_cst::string_type>
void execute(const char* input, uint8_t num_bytes, t_pat pat, const char* format);

    //const char* input = "2 801 543 293 597 801 444 444 293";
    const char* input = "1 2 1 4 3 4 1 4 2 2 2 1 1 1 4 2 4";
    uint8_t num_bytes = 'd';
    std::cout << (int) num_bytes << std::endl;
    execute<CST_Type>(input, num_bytes, {801, 444}, "%2I %3S %:4T");

    std::cout << "DFS: " << std::endl;
    CST_Type cst;
    sdsl::construct_im(cst, input, num_bytes);
    traverse(cst);

}


void traverse(CST_Type &cst) {
    uint64_t max_depth = 200;

    // use the DFS iterator to traverse `cst`
    for (auto it=cst.begin(); it!=cst.end(); ++it) {
        if (it.visit() == 1) {  // node visited the first time
            auto v = *it;       // get the node by dereferencing the iterator
            if (cst.depth(v) <= max_depth) {   // if depth node is <= max_depth
                // process node, e.g. output it in format d-[lb, rb]
                std::cout<<cst.depth(v)<<"-["<<cst.lb(v)<< ","<<cst.rb(v)<<"]"<<std::endl;
            } else { // skip the subtree otherwise
                it.skip_subtree();
            }
        }
    }
}
*/


std::string cout_vector(const std::vector<long unsigned int> &v) {
    std::stringstream s;
    for (const auto &i: v)
        s << i << ", ";
    return s.str();
}


template<class t_cst, class t_pat=typename t_cst::string_type>
void execute(const char *input, uint8_t num_bytes, t_pat pat, const char *format) {
    typedef typename t_cst::node_type node_t;
    t_cst cst;
    sdsl::construct_im(cst, input, num_bytes);
    sdsl::csXprintf(std::cout, format, cst);

    std::cout << "pattern \"" << cout_vector(pat) << "\"" << std::endl;
    std::cout << "---- backward search step by step ----" << std::endl;
    {
        uint64_t lb = 0, rb = cst.size() - 1;
        for (auto it = pat.end(); it != pat.begin() and lb <= rb;) {
            --it;
            if (sdsl::backward_search(cst.csa, lb, rb, (typename t_cst::char_type) *it, lb, rb) > 0) {
                std::cout << "[" << lb << "," << rb << "]" << std::endl;
                std::cout << "matched " << *it << std::endl;
            }
        }
    }

    std::cout << "---- backward search for the pattern  ----" << std::endl;
    {
        uint64_t lb = 0, rb = cst.size() - 1;
        sdsl::backward_search(cst.csa, lb, rb, pat.begin(), pat.end(), lb, rb);
        std::cout << "size = " << rb + 1 - lb << std::endl;
    }
    std::cout << "---- count pattern occurrences  ----" << std::endl;
    {
        std::cout << "count(cst.csa, \"" << cout_vector(pat) << "\")=" << sdsl::count(cst.csa, pat) << std::endl;
    }
    std::cout << "---- locate the pattern  ----" << std::endl;
    {
        auto occs = sdsl::locate(cst.csa, pat);
        std::cout << "locate(cst.csa, \"" << cout_vector(pat) << "\")=" << occs.size() << std::endl;
        std::cout << occs << std::endl;
        std::cout << std::endl;
    }
    std::cout << "---- extract text  ----" << std::endl;
    {
        std::cout << "extract(csa,0,csa.size())=\"" << cout_vector(sdsl::extract(cst.csa, 0, cst.csa.size() - 1))
                  << "\"" << std::endl;
    }

    std::cout << "---- forward search step by step  ----" << std::endl;
    {
        node_t v = cst.root();
        auto it = pat.begin();
        for (uint64_t char_pos = 0; it != pat.end(); ++it) {
            if (sdsl::forward_search(cst, v, it - pat.begin(), *it, char_pos) > 0) {
                std::cout << it - pat.begin() << "-[" << cst.lb(v) << "," << cst.rb(v) << "]" << std::endl;
                std::cout << "matched " << *it << std::endl;
            } else {
                break;
            }
        }
    }
    std::cout << "---- count pattern occurrences ----" << std::endl;
    {
        std::cout << "count(cst, \"" << cout_vector(pat) << "\")=" << sdsl::count(cst, pat) << std::endl;
    }
    std::cout << "---- extract text  ----" << std::endl;
    {
        std::cout << "extract(cst,cst.select_leaf(cst.csa.isa[0]+1))=\""
                  << cout_vector(sdsl::extract(cst, cst.select_leaf(cst.csa.isa[0] + 1))) << "\"" << std::endl;
    }

}