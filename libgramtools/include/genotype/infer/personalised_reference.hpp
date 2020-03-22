#ifndef GRAM_PERSONALISED_REF_H
#define GRAM_PERSONALISED_REF_H

#include <string>
#include <vector>
#include "prg/types.hpp"
#include "genotype/infer/types.hpp"

using namespace gram::genotype::infer;

namespace gram::genotype{
    class SegmentTracker;

    static constexpr int FASTA_LWIDTH = 60;
    class Fasta{
        std::string ID{""}, desc{""};
        std::string sequence;
    public:
        std::string const& get_sequence() {return sequence;}
        void set_ID(std::string new_ID) {this->ID = new_ID;}
        void set_desc(std::string new_desc) {this->desc = new_desc;}
        void add_sequence(std::string const& seq) {sequence += seq;}

        friend bool operator<(const Fasta& first, const Fasta& second);
        friend std::ostream& operator<<(std::ostream& out_stream, const Fasta& input);
    };

    using Fastas = std::vector<Fasta>;
    using unique_Fastas = std::set<Fasta>;


    class InconsistentPloidyException : public std::exception{
    public:
        InconsistentPloidyException() {}
        virtual char const* what() const throw(){
            return "The sites do not all have the same GT cardinality (ploidy)"; }
    };

    allele_vector get_all_alleles_to_paste(gt_site_ptr const& site, std::size_t ploidy);
    Fastas get_personalised_ref(covG_ptr graph_root, gt_sites const &genotyped_records,
            SegmentTracker &tracker);

    void add_description(Fastas &p_refs, std::string const &desc);
}

#endif //GRAM_PERSONALISED_REF_H
