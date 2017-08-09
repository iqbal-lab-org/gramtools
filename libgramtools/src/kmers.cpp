#include <algorithm>
#include <thread>

#include "fm_index.hpp"
#include "bidir_search_bwd.hpp"
#include "kmers.hpp"

#define MAX_THREADS 25


inline bool fexists(const std::string &name) {
    std::ifstream f(name.c_str());
    return f.good();
}


//trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}


// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}


// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}


std::vector<std::string> split(const std::string &cad, const std::string &delim) {
    std::vector<std::string> v;
    uint64_t p, q, d;

    q = cad.size();
    p = 0;
    d = strlen(delim.c_str());
    while (p < q) {
        uint64_t posfound = cad.find(delim, p);
        std::string token;

        if (posfound >= 0)
            token = cad.substr(p, posfound - p);
        else
            token = cad.substr(p, cad.size() + 1);
        p += token.size() + d;
        trim(token);
        v.push_back(token);
    }
    return v;
}


void calc_kmer_matches(KmerIdx &kmer_idx, KmerIdx &kmer_idx_rev,
                       KmerSites &kmer_sites,
                       sequence_set<std::vector<uint8_t>> &kmers_in_ref,
                       std::vector<std::vector<uint8_t>> &kmers,
                       const uint64_t maxx,
                       const std::vector<int> &mask_a,
                       const DNA_Rank &rank_all,
                       const FM_Index &fm_index,
                       const int thread_id) {

    for (auto &kmer: kmers) {
        kmer_idx[kmer] = std::list<std::pair<uint64_t, uint64_t>>();
        kmer_idx_rev[kmer] = std::list<std::pair<uint64_t, uint64_t>>();
        kmer_sites[kmer] = std::list<std::vector<std::pair<uint64_t, std::vector<int>>>>();

        bool delete_first_interval = false;
        bool kmer_precalc_done = false;

        bidir_search_bwd(kmer_idx[kmer], kmer_idx_rev[kmer],
                         0, fm_index.size(),
                         0, fm_index.size(),
                         kmer_sites[kmer],
                         delete_first_interval,
                         kmer.begin(), kmer.end(),
                         mask_a, maxx, kmer_precalc_done,
                         rank_all, fm_index, thread_id);

        if (kmer_idx[kmer].empty())
            kmer_idx.erase(kmer);

        if (kmer_idx_rev[kmer].empty())
            kmer_idx_rev.erase(kmer);

        if (!delete_first_interval)
            kmers_in_ref.insert(kmer);
    }
}


void *worker(void *st) {
    auto *th = (ThreadData *) st;
    calc_kmer_matches(*(th->kmer_idx),
                      *(th->kmer_idx_rev),
                      *(th->kmer_sites),
                      *(th->kmers_in_ref),
                      *(th->kmers),
                      th->maxx,
                      *(th->mask_a),
                      *(th->rank_all),
                      *(th->fm_index),
                      th->thread_id);
    return nullptr;
}


void generate_kmers_encoding(const std::vector<int> &mask_a,
                             const std::string &kmer_fname,
                             const uint64_t maxx, const int k,
                             const int thread_count,
                             const DNA_Rank &rank_all,
                             const FM_Index &fm_index) {

    pthread_t threads[thread_count];
    struct ThreadData td[thread_count];
    std::vector<std::vector<uint8_t>> kmers[thread_count];

    std::ifstream kfile;
    std::string line;
    kfile.open(kmer_fname);

    int j = 0;
    while (std::getline(kfile, line)) {
        std::vector<uint8_t> kmer;
        for (const auto &c: line) {
            switch (c) {
                case 'A':
                case 'a':
                    kmer.push_back(1);
                    break;
                case 'C':
                case 'c':
                    kmer.push_back(2);
                    break;
                case 'G':
                case 'g':
                    kmer.push_back(3);
                    break;
                case 'T':
                case 't':
                    kmer.push_back(4);
                    break;
                default:
                    break;
            }
        }
        kmers[j++].push_back(kmer);
        j %= thread_count;
    }

    sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>> kmer_idx[thread_count], kmer_idx_rev[thread_count];
    sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint64_t, std::vector<int>>>>> kmer_sites[thread_count];
    sequence_set<std::vector<uint8_t>> kmers_in_ref[thread_count];

    for (int i = 0; i < thread_count; i++) {
        td[i].fm_index = &fm_index;
        td[i].k = k;
        td[i].kmer_idx = &kmer_idx[i];
        td[i].kmer_idx_rev = &kmer_idx_rev[i];
        td[i].kmer_sites = &kmer_sites[i];
        td[i].mask_a = &mask_a;
        td[i].maxx = maxx;
        td[i].kmers_in_ref = &kmers_in_ref[i];
        td[i].kmers = &kmers[i];
        td[i].thread_id = i;
        std::cout << "Starting thread: " << i << std::endl;
        pthread_create(&threads[i], nullptr, worker, &td[i]);
    }

    std::ofstream precalc_file;
    precalc_file.open(std::string(kmer_fname) + ".precalc");

    std::cout << "Joining threads" << std::endl;

    for (int i = 0; i < thread_count; i++) {
        void *status;
        pthread_join(threads[i], &status);
        for (auto obj: kmer_idx[i]) {
            auto k = obj.first;
            for (const auto &n: k)
                precalc_file << (int) n << ' ';

            precalc_file << '|';

            if (kmers_in_ref[i].count(k))
                precalc_file << 1;
            else
                precalc_file << 0;

            precalc_file << '|';
            for (auto o: kmer_idx[i][k])
                precalc_file << o.first << ' ' << o.second << ' ';
            precalc_file << '|';
            for (auto o: kmer_idx_rev[i][k])
                precalc_file << o.first << ' ' << o.second << ' ';
            precalc_file << '|';

            for (const auto &o: kmer_sites[i][k]) {
                for (const auto &v: o) {
                    precalc_file << v.first << ' ';
                    for (const auto &n: v.second)
                        precalc_file << (int) n << ' ';
                    precalc_file << '@';
                }
                precalc_file << '|';
            }

            precalc_file << std::endl;
        }
    }
}


KmersData read_encoded_kmers(const std::string &encoded_kmers_fname) {
    std::ifstream fhandle;
    fhandle.open(encoded_kmers_fname);

    KmersData kmers;
    std::string line;
    while (std::getline(fhandle, line)) {

        std::vector<std::string> parts = split(line, "|");
        std::vector<uint8_t> kmer;
        for (const auto &c: split(parts[0], " "))
            kmer.push_back(std::stoi(c));

        if (parts[1] != "1")
            kmers.in_reference.insert(kmer);

        std::vector<std::string> idx_str = split(parts[2], " ");
        std::vector<std::string> idx_rev_str = split(parts[3], " ");
        std::list<std::pair<uint64_t, uint64_t>> idx, idx_rev;
        for (unsigned int i = 0; i < idx_str.size() / 2; i++)
            idx.emplace_back(std::pair<uint64_t, uint64_t>(std::stoi(idx_str[i * 2]),
                                                           std::stoi(idx_str[i * 2 + 1])));
        for (unsigned int i = 0; i < idx_rev_str.size() / 2; i++)
            idx_rev.emplace_back(
                    std::pair<uint64_t, uint64_t>(std::stoi(idx_rev_str[i * 2]),
                                                  std::stoi(idx_rev_str[i * 2 + 1])));

        int flag = 0;
        if (!idx.empty()) {
            kmers.index[kmer] = idx;
            flag = 1;
        }
        if (!idx_rev.empty())
            kmers.index_reverse[kmer] = idx_rev;

        if (flag != 1)
            continue;

        std::list<std::vector<std::pair<uint64_t, std::vector<int>>>> sites;
        for (unsigned int i = 4; i < parts.size(); i++) {
            std::vector<std::pair<uint64_t, std::vector<int>>> site;
            for (auto pair_i_v: split(parts[i], "@")) {
                std::vector<std::string> strvec = split(pair_i_v, " ");
                if (strvec.size()) {
                    int first = std::stoi(strvec[0]);

                    std::vector<int> rest;
                    for (unsigned int j = 1; j < strvec.size(); j++)
                        if (strvec[j].size())
                            rest.push_back(std::stoi(strvec[j]));

                    site.emplace_back(std::pair<uint64_t, std::vector<int>>(first, rest));
                }
            }
            sites.push_back(site);
        }
        kmers.sites[kmer] = sites;
    }
    return kmers;
}


int get_thread_count() {
    int count_total = std::thread::hardware_concurrency();
    int selected_count = std::min(MAX_THREADS, count_total) - 1;
    if (selected_count < 1)
        return 1;
    return selected_count;
}


KmersData get_kmers(const std::string &kmer_fname,
                    const int k,
                    const std::vector<int> &mask_a,
                    const uint64_t maxx,
                    const DNA_Rank &rank_all,
                    const FM_Index &fm_index) {

    const auto encoded_kmers_fname = std::string(kmer_fname) + ".precalc";

    if (!fexists(encoded_kmers_fname)) {
        const int thread_count = get_thread_count();
        std::cout << "Precalculated kmers not found, generating them using "
                  << thread_count << " threads" << std::endl;
        generate_kmers_encoding(mask_a, kmer_fname, maxx, k, thread_count, rank_all, fm_index);
        std::cout << "Finished precalculating kmers" << std::endl;
    }

    std::cout << "Reading Kmers" << std::endl;
    const auto kmers = read_encoded_kmers(encoded_kmers_fname);
    return kmers;
}
