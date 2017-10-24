#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "sequence_read/seqread.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "search.hpp"
#include "coverage_analysis.hpp"


uint64_t quasimap_reads(const Parameters &params,
                        const KmerIndex &kmer_index,
                        const PRG_Info &prg_info) {
    auto allele_coverage = make_allele_coverage_structure(prg_info);

    SeqRead reads(params.reads_fpath.c_str());
    std::ofstream progress_fhandle(params.processed_reads_fpath);

    uint64_t count_all_reads = 0;
    uint64_t count_mapped_reads = 0;

    for (const auto *const raw_read: reads) {
        if (count_all_reads++ % 100000 == 0)
            progress_fhandle << count_all_reads << std::endl;

        auto read = encode_dna_bases(*raw_read);
        bool read_mappred_exactly = quasimap_read(read,
                                                  allele_coverage,
                                                  kmer_index,
                                                  prg_info,
                                                  params);
        if (read_mappred_exactly)
            ++count_mapped_reads;
    }
    progress_fhandle.close();
    return count_mapped_reads;
}


bool quasimap_read(const Pattern &read,
                   AlleleCoverage &allele_coverage,
                   const KmerIndex &kmer_index,
                   const PRG_Info &prg_info,
                   const Parameters &params) {
    const auto kmer = get_kmer_from_read(params.kmers_size, read);
    const auto search_states = search_read_bwd(read, kmer, kmer_index, prg_info);
    const bool read_mapped_exactly = not search_states.empty();

    record_read_coverage(allele_coverage, search_states);
    return read_mapped_exactly;
}


void record_read_coverage(AlleleCoverage &allele_coverage,
                          const SearchStates &search_states) {
    for (const auto &search_state: search_states) {
        for (const auto &variant_site: search_state.variant_site_path) {
            auto marker = variant_site.first;
            auto allell_id = variant_site.second;

            auto min_boundary_marker = 5;
            auto variant_site_cover_index = (marker - min_boundary_marker) / 2;
            auto allele_cover_index = allell_id - 1;

            allele_coverage[variant_site_cover_index][allele_cover_index] += 1;
        }
    }
}


AlleleCoverage make_allele_coverage_structure(const PRG_Info &prg_info) {
    AlleleCoverage allele_coverage;

    Marker last_variant_site_marker = 0;
    uint32_t count_sites = 0;

    for (const auto &prg_char: prg_info.fm_index.text) {
        if (prg_char <= 4)
            continue;

        bool char_is_allele_boundary_marker = prg_char % 2 == 0;
        if (char_is_allele_boundary_marker) {
            count_sites++;
            continue;
        }

        auto done_with_variant_site = last_variant_site_marker == prg_char;
        if (done_with_variant_site) {
            std::vector<uint32_t> allele(count_sites + 1);
            std::fill(allele.begin(), allele.end(), 0);
            allele_coverage.emplace_back(allele);
            count_sites = 0;
            continue;
        }

        last_variant_site_marker = prg_char;

    }
    return allele_coverage;
}


Pattern get_kmer_from_read(const uint32_t kmer_size, const Pattern &read) {
    Pattern kmer;
    auto kmer_start_it = read.begin() + read.size() - kmer_size;
    kmer.assign(kmer_start_it, read.end());
    return kmer;
}
