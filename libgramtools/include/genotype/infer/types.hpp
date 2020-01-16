#ifndef GT_INFER_TYPES
#define GT_INFER_TYPES

#include <memory>
#include "prg/coverage_graph.hpp"

namespace gram::genotype::infer {
    class AbstractGenotypedSite;
    using gt_site = AbstractGenotypedSite;
    using gt_site_ptr = std::shared_ptr<AbstractGenotypedSite>;
    using gt_sites = std::vector<gt_site_ptr>;

    struct Allele {
        std::string sequence;
        PerBaseCoverage pbCov;
        AlleleId haplogroup; /**< Which ID in its site this allele is associated with */

        /**
         * Allele combination overload
         * The left-hand side (= this object) argument's haplogroup is used, regardless of `other`'s haplogroup.
         */
        Allele operator+(Allele const &other) const {
            PerBaseCoverage new_pbCov;
            new_pbCov.reserve(pbCov.size() + other.pbCov.size());
            new_pbCov.insert(new_pbCov.end(), pbCov.begin(), pbCov.end());
            new_pbCov.insert(new_pbCov.end(), other.pbCov.begin(), other.pbCov.end());
            return Allele{
                    sequence + other.sequence,
                    new_pbCov,
                    haplogroup
            };
        }

        bool operator==(Allele const& other) const{
            return sequence == other.sequence &&
                   pbCov == other.pbCov &&
                   haplogroup == other.haplogroup;
        }
    };

    using allele_vector = std::vector<Allele>;

    enum class Ploidy{Haploid, Diploid};

}

#endif //GT_INFER_TYPES
