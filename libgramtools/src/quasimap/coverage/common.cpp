#include "quasimap/coverage/allele_sum.hpp"
#include "quasimap/coverage/allele_base.hpp"
#include "quasimap/coverage/grouped_allele_counts.hpp"

#include "quasimap/coverage/common.hpp"


void coverage::record::search_states(Coverage &coverage,
                                     const SearchStates &search_states,
                                     const uint64_t &read_length,
                                     const PRG_Info &prg_info) {
    coverage::record::allele_sum(coverage, search_states);
    coverage::record::grouped_allele_counts(coverage, search_states);
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);
}


void coverage::dump(const Coverage &coverage,
                   const Parameters &parameters) {
    std::ofstream file_handle(parameters.allele_coverage_fpath);
    for (const auto &variant_site_coverage: coverage.allele_sum_coverage) {
        auto allele_count = 0;
        for (const auto &sum_coverage: variant_site_coverage) {
            file_handle << sum_coverage;
            auto not_last_coverage = allele_count++ < variant_site_coverage.size() - 1;
            if (not_last_coverage)
                file_handle << " ";
        }
        file_handle << std::endl;
    }
}


Coverage coverage::generate::empty_structure(const PRG_Info &prg_info) {
    Coverage coverage;
    coverage.allele_sum_coverage = coverage::generate::allele_sum_structure(prg_info);
    coverage.allele_base_coverage = coverage::generate::allele_base_structure(prg_info);
    coverage.grouped_allele_counts = coverage::generate::grouped_allele_counts(prg_info);
    return coverage;
}