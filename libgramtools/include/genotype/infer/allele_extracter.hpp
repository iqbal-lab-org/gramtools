#include "common/utils.hpp"
#include "types.hpp"
#include "prg/coverage_graph.hpp"

using namespace gram;

namespace gram::infer{

struct Allele{
    std::string sequence;
    BaseCoverage pbCov;

    int haplogroup; /**< Which ID in its site this allele is associated with */

    Allele operator+(Allele const& other) const{
        BaseCoverage new_pbCov(pbCov.size() + other.pbCov.size());
        new_pbCov.insert(new_pbCov.end(), pbCov.begin(), pbCov.end());
        new_pbCov.insert(new_pbCov.end(), other.pbCov.begin(), other.pbCov.end());
        return Allele{
           sequence + other.sequence,
           new_pbCov
        };
    }

    void append(std::string const& additional_seq, BaseCoverage const& additional_pbCov){
       sequence += additional_seq;
       pbCov.insert( pbCov.end(), additional_pbCov.begin(), additional_pbCov.end() );
    }
};


/**
 * Class in charge of producing the set of `Allele`s that get genotyped.
 * The procedure scans through each haplogroup of a site, pasting sequence & coverage
 * from previously genotyped (=nested) sites when encountered.
 */
class AlleleExtracter{
private:
    allele_vector alleles;
    gt_sites* genotyped_sites;
public:
    AlleleExtracter() = default;
    AlleleExtracter(covG_ptr site_start, covG_ptr site_end, gt_sites const& genotyped_sites);

    allele_vector get_alleles() const{return alleles;}

    /**
     * Linear traversal of an allelic haplogroup, extracting all relevant combinations of alleles.
     * In the absence of incident nested sites, always produces a single allele.
     */
    allele_vector extract_alleles(covG_ptr haplogroup_start, covG_ptr site_end);

    /**
     * From the set of genotypes of a site, makes new alleles containing the genotyped alleles.
     * @param existing the starting alleles
     * @param site_index the previously genotyped site
     * @return Cartesian product of `existing` and distinct genotyped alleles
     */
    allele_vector allele_combine(allele_vector existing, std::size_t site_index);
};

