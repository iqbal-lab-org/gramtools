#ifndef GT_INFER_TYPES
#define GT_INFER_TYPES

#include <memory>

namespace gram::genotype::infer {
    class AbstractGenotypedSite;
    using gtSite_ptr = std::shared_ptr<AbstractGenotypedSite>;

    struct Allele;
    using gt_site = AbstractGenotypedSite;
    using gt_sites = std::vector<gtSite_ptr>;
    using allele_vector = std::vector<Allele>;
}

#endif //GT_INFER_TYPES
