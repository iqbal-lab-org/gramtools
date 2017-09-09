#include <list>

#include "bidir_search_bwd.hpp"
#include "kmer_index.hpp"


using Kmer = std::vector<uint8_t>;
using Kmers = std::vector<Kmer>;

using SA_Interval = std::pair<uint64_t, uint64_t>;
using SA_Intervals = std::list<SA_Interval>;


void print_sa_intervals(const SA_Intervals &sa_intervals) {
    for (const auto &sa_interval: sa_intervals)
        std::cout << "(" << sa_interval.first << ", " << sa_interval.second << ")" << "   ";
    std::cout << std::endl;
}


struct CacheElement {
    SA_Intervals sa_intervals;
    Sites sites;
    uint8_t base = 0;
};


using KmerCache = std::list<CacheElement>;


void print_cache(const KmerCache &cache) {
    for (const auto element: cache) {
        std::cout << (int) element.base << ", ";
    }
    std::cout << std::endl;
}


void update_cache(std::list<CacheElement> &cache,
                  const Kmer &kmer_suffix_diff,
                  const uint64_t kmer_size,
                  const uint64_t max_alphabet_num,
                  const std::vector<int> &allele_mask,
                  const DNA_Rank &rank_all,
                  const FM_Index &fm_index) {

    auto new_size = kmer_size - kmer_suffix_diff.size();
    cache.resize(new_size);

    for (const auto &base: kmer_suffix_diff) {
        CacheElement new_cache_element;
        new_cache_element.base = base;

        if (cache.size() == 0) {
            new_cache_element.sa_intervals = {{0, fm_index.size()}};
            new_cache_element.sites = {Site()};
        } else {
            const auto &last_cache_element = cache.back();
            new_cache_element.sa_intervals = last_cache_element.sa_intervals;
            new_cache_element.sites = last_cache_element.sites;
        }

        reduce_sa_intervals(base,
                            new_cache_element.sa_intervals,
                            new_cache_element.sites,
                            false, true,
                            allele_mask,
                            max_alphabet_num,
                            false,
                            rank_all,
                            fm_index);

        cache.emplace_back(new_cache_element);
    }
}


void generate_kmer_index(const Kmers &kmer_suffix_diffs,
                         const uint64_t max_alphabet_num,
                         const std::vector<int> &allele_mask,
                         const DNA_Rank &rank_all,
                         const FM_Index &fm_index) {

    const uint64_t kmer_size = kmer_suffix_diffs[0].size();
    std::list<CacheElement> cache;

    for (const auto &kmer_suffix_diff: kmer_suffix_diffs) {
        update_cache(cache,
                     kmer_suffix_diff,
                     kmer_size,
                     max_alphabet_num,
                     allele_mask,
                     rank_all,
                     fm_index);

        print_cache(cache);

        for (const auto &elem: cache)
            print_sa_intervals(elem.sa_intervals);

        //break;

        // do things with last cache element
    }
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