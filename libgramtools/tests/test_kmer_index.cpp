#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <sdsl/suffix_arrays.hpp>

#include "gtest/gtest.h"

#include "utils.hpp"
#include "map.hpp"
#include "prg.hpp"
#include "kmer_index.hpp"


class GenerateKmerIndex : public ::testing::Test {

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


TEST_F(GenerateKmerIndex, hip) {
    const std::string prg_raw = "catttacatt"
            "5c6t5"
            "aaagcaacagaac";
    const uint64_t max_alphabet = max_alphabet_num(prg_raw);
    const std::vector<int> allele_mask = generate_allele_mask(prg_raw);

    const FM_Index fm_index = fm_index_from_raw_prg(prg_raw);
    const DNA_Rank &rank_all = calculate_ranks(fm_index);

    Kmers kmers = {
            {1, 2, 2},
            {1, 2, 3},
            {1, 3, 1},
            /*
            encode_dna_bases("catt"),
            encode_dna_bases("attt"),
            encode_dna_bases("ttta"),
             */
    };

    generate_kmer_index(kmers, max_alphabet, allele_mask, rank_all, fm_index);
}
