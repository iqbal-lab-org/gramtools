#ifndef GT_INFER_TYPES
#define GT_INFER_TYPES

#include <memory>

#include "common/data_types.hpp"

namespace gram {
struct Allele {
  std::string sequence;
  PerBaseCoverage pbCov;
  AlleleId
      haplogroup; /**< Which ID in its site this allele is associated with */
  bool callable = true; /**< Whether the allele is consistent with
                                       calls in child sites */

  Allele() : sequence(""), haplogroup(0) {}
  Allele(std::string seq, PerBaseCoverage pb)
      : sequence(seq), pbCov(pb), haplogroup(0) {}
  Allele(std::string seq, PerBaseCoverage pb, AlleleId haplo)
      : sequence(seq), pbCov(pb), haplogroup(haplo) {}
  Allele(std::string seq, PerBaseCoverage pb, AlleleId haplo, bool consistent)
      : sequence(seq), pbCov(pb), haplogroup(haplo), callable(consistent) {}

  /**
   * Allele combination overload
   * The left-hand side (= this object) argument's haplogroup is used,
   * regardless of `other`'s haplogroup.
   * The right-hand site `callable` if false always overrides,because
   * any non-callable portion makes whole allele uncallable.
   */
  Allele operator+(Allele const& other) const {
    PerBaseCoverage new_pbCov;
    new_pbCov.reserve(pbCov.size() + other.pbCov.size());
    new_pbCov.insert(new_pbCov.end(), pbCov.begin(), pbCov.end());
    new_pbCov.insert(new_pbCov.end(), other.pbCov.begin(), other.pbCov.end());
    bool new_consistent = callable & other.callable;
    return Allele{sequence + other.sequence, new_pbCov, haplogroup,
                  new_consistent};
  }

  bool operator<(Allele const& other) const {
    return sequence < other.sequence;
  }

  /**
   * Note: callable boolean attribute should NOT be used
   * in comparison, e.g., when looking for the ref in vector of
   * alleles, either state should give a match.
   */
  bool operator==(Allele const& other) const {
    return sequence == other.sequence && pbCov == other.pbCov &&
           haplogroup == other.haplogroup;
  }

  double get_average_cov() const {
    double result = std::accumulate(pbCov.begin(), pbCov.end(), 0.0);
    result /= pbCov.size();
    return result;
  }
};

using allele_vector = std::vector<Allele>;
}  // namespace gram

namespace gram::genotype::infer {
class Genotyper;
using gtyper_ptr = std::shared_ptr<Genotyper>;

class GenotypedSite;
using gt_site = GenotypedSite;
using gt_site_ptr = std::shared_ptr<GenotypedSite>;
using gt_sites = std::vector<gt_site_ptr>;

using GtypedIndex =
    int32_t; /**< The index of an allele in an allele vector or -1 for null*/
using GtypedIndices = std::vector<GtypedIndex>;
using allele_coverages = std::vector<double>;

}  // namespace gram::genotype::infer

#endif  // GT_INFER_TYPES
