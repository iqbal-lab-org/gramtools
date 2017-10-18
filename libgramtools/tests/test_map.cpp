#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <sdsl/suffix_arrays.hpp>

#include "gtest/gtest.h"
#include "map.hpp"


class QuasimapRead : public ::testing::Test {

protected:
    std::string prg_fpath;

    virtual void SetUp() {
        boost::uuids::uuid uuid = boost::uuids::random_generator()();
        const auto uuid_str = boost::lexical_cast<std::string>(uuid);
        prg_fpath = "./prg_" + uuid_str;
    }

    virtual void TearDown() {
        std::remove(prg_fpath.c_str());
    }

    FM_Index fm_index_from_raw_prg(const std::string &prg_raw) {
        std::vector<uint64_t> prg = encode_prg(prg_raw);
        dump_encoded_prg(prg, prg_fpath);
        FM_Index fm_index;
        // TODO: constructing from memory with sdsl::construct_im appends 0 which corrupts
        sdsl::construct(fm_index, prg_fpath, 8);
        return fm_index;
    }
};

/*
TEST_F(QuasimapRead, todo_desc) {
    const std::string prg_raw = "taca5g6t5aat";
    const std::string read = "acagaat";
    const std::vector<uint8_t> kmer = encode_read("agaat");

    Patterns kmers = {kmer};

    MasksParser masks;
    masks.allele = generate_allele_mask(prg_raw);
    masks.max_alphabet_num = max_alphabet_num(prg_raw);

    const auto fm_index = fm_index_from_raw_prg(prg_raw);
    const auto rank_all = calculate_ranks(fm_index);
    const auto encoded_read = encode_dna_bases(read);

    KmerIndex kmers_data;
    index_kmers(kmers,
                kmers_data.sa_intervals_map,
                kmers_data.sites_map,
                kmers_data.non_site_crossing_kmers,
                masks.max_alphabet_num,
                masks.allele,
                rank_all,
                fm_index);

    for (auto i: kmers_data.sites_map) {
        std::cout << "kmer: " << std::endl;
        auto &kmer = i.first;
        for (auto base: kmer)
            std::cout << (int) base << " ";
        std::cout << std::endl;

        std::cout << "site: " << std::endl;
        auto &sites = i.second;
        print_sites(sites);

        std::cout << "sa_intervals: " << std::endl;
        auto &sa_intervals = kmers_data.sa_intervals_map[kmer];
        print_sa_interval(sa_intervals);
    }

    int count_char_in_variant_site = 0;
    std::unordered_set<uint64_t> repeats_variant_site_edge_markers;

    bool read_mapped = quasimap_read(kmer,
                                     encoded_read,
                                     count_char_in_variant_site,
                                     repeats_variant_site_edge_markers,
                                     kmers_data,
                                     masks,
                                     kmer.size(),
                                     rank_all,
                                     fm_index);

    std::cout << "size: " <<  masks.allele_coverage.size() << std::endl;

    for (uint64_t i = 0; i < masks.allele_coverage.size(); i++) {
        for (uint64_t j = 0; j < masks.allele_coverage[i].size(); j++)
            std::cout << masks.allele_coverage[i][j] << " ";
        std::cout << std::endl;
    }
}
*/