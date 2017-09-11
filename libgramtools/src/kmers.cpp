#include <algorithm>
#include <thread>

#include "utils.hpp"
#include "fm_index.hpp"
#include "kmers.hpp"
#include "bidir_search_bwd.hpp"


std::string dump_kmer(const Kmer &kmer) {
    std::stringstream stream;
    for (const auto &base: kmer) {
        stream << (int) base;
        bool last_base = (&base == &kmer.back());
        if (!last_base)
            stream << ' ';
    }
    return stream.str();
}


std::string dump_sa_intervals(const SA_Intervals &sa_intervals) {
    std::stringstream stream;
    for (const auto &sa_interval: sa_intervals) {
        stream << sa_interval.first
               << " "
               << sa_interval.second;
        bool last_base = (&sa_interval == &sa_intervals.back());
        if (!last_base)
            stream << " ";
    }
    return stream.str();
}


std::string dump_crosses_marker_flag(const Kmer &kmer,
                                     const NonVariantKmers &nonvar_kmers) {
    if (nonvar_kmers.count(kmer) != 0)
        return "1";
    return "0";
}


std::string dump_sites(const Kmer &kmer, const KmerSites &kmer_sites) {
    std::stringstream stream;
    for (const auto &all_sites: kmer_sites.at(kmer)) {
        for (const auto &sites: all_sites) {
            stream << sites.first << " ";
            for (const auto &allele: sites.second)
                stream << (int) allele << " ";
            stream << "@";
        }
        stream << "|";
    }
    return stream.str();
}


std::string dump_kmer_index_entry(const Kmer &kmer,
                                  const SA_Intervals &sa_intervals,
                                  const NonVariantKmers &nonvar_kmers,
                                  const KmerSites &kmer_sites) {
    std::stringstream stream;
    stream << dump_kmer(kmer) << "|";
    stream << dump_crosses_marker_flag(kmer, nonvar_kmers) << "|";
    // additional bar because of reverse sa intervals
    stream << dump_sa_intervals(sa_intervals) << "||";
    stream << dump_sites(kmer, kmer_sites);
    return stream.str();
}


void dump_kmer_index(std::ofstream &kmer_index_file,
                     const KmerSA_Intervals &kmers_sa_intervals,
                     const NonVariantKmers &nonvar_kmers,
                     const KmerSites &kmer_sites) {

    for (const auto &kmer_sa_intervals: kmers_sa_intervals) {
        const auto &kmer = kmer_sa_intervals.first;
        const SA_Intervals &sa_intervals = kmers_sa_intervals.at(kmer);
        const std::string kmer_entry =
                dump_kmer_index_entry(kmer,
                                      sa_intervals,
                                      nonvar_kmers,
                                      kmer_sites);
        kmer_index_file << kmer_entry << std::endl;
    }
}


void index_kmers(Kmers &kmers,
                 KmerSA_Intervals &sa_intervals_map,
                 KmerSites &sites_map,
                 NonVariantKmers &nonvar_kmers,
                 const uint64_t max_alphabet_num,
                 const std::vector<int> &allele_mask,
                 const DNA_Rank &rank_all,
                 const FM_Index &fm_index) {

    for (auto &kmer: kmers) {
        auto &sa_intervals = sa_intervals_map[kmer];
        auto &sites = sites_map[kmer];

        sa_intervals = SA_Intervals();
        sites = Sites();

        bool delete_first_interval = false;
        bool kmer_precalc_done = false;

        if (sa_intervals.empty()) {
            sa_intervals.emplace_back(std::make_pair(0, fm_index.size()));
            sites.emplace_back(Site());
        }

        bidir_search_bwd(sa_intervals, sites,
                         delete_first_interval,
                         kmer.begin(),
                         kmer.end(),
                         allele_mask,
                         max_alphabet_num,
                         kmer_precalc_done,
                         rank_all, fm_index);

        if (sa_intervals.empty()) {
            sa_intervals_map.erase(kmer);
        }

        if (!delete_first_interval)
            nonvar_kmers.insert(kmer);
    }
}


void generate_kmer_index(const std::vector<int> &allele_mask,
                         const std::string &kmer_fname,
                         const uint64_t max_alphabet_num,
                         const DNA_Rank &rank_all,
                         const FM_Index &fm_index) {

    std::ifstream kmer_fhandle;
    kmer_fhandle.open(kmer_fname);

    Kmers kmers;
    std::string line;
    while (std::getline(kmer_fhandle, line)) {
        const Kmer &kmer = encode_dna_bases(line);
        kmers.emplace_back(kmer);
    }

    KmerSA_Intervals sa_intervals_map;
    KmerSites sites_map;
    NonVariantKmers nonvar_kmers;

    index_kmers(kmers,
                sa_intervals_map, sites_map, nonvar_kmers,
                max_alphabet_num, allele_mask, rank_all,
                fm_index);

    std::ofstream kmer_index_file;
    kmer_index_file.open(std::string(kmer_fname) + ".precalc");

    dump_kmer_index(kmer_index_file,
                    sa_intervals_map,
                    nonvar_kmers,
                    sites_map);
}


inline bool file_exists(const std::string &name) {
    std::ifstream f(name.c_str());
    return f.good();
}


static inline std::string &left_trim(std::string &s) {
    //trim from start
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}


static inline std::string &right_trim(std::string &s) {
    // trim from end
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}


static inline std::string &trim(std::string &s) {
    // trim from both ends
    return left_trim(right_trim(s));
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


bool parse_crosses_marker_flag(const std::string &in_reference_flag_str) {
    return in_reference_flag_str == "1";
}


Kmer parse_encoded_kmer(const std::string &encoded_kmer_str) {
    Kmer encoded_kmer;
    for (const auto &encoded_base: split(encoded_kmer_str, " ")) {
        const auto kmer_base = (uint8_t) std::stoi(encoded_base);
        encoded_kmer.push_back(kmer_base);
    }
    return encoded_kmer;
}


SA_Intervals parse_sa_intervals(const std::string &full_sa_intervals_str) {
    std::vector<std::string> split_sa_intervals = split(full_sa_intervals_str, " ");
    SA_Intervals sa_intervals;
    for (auto i = 0; i < split_sa_intervals.size(); i = i + 2) {
        auto sa_interval_start = (uint64_t) std::stoi(split_sa_intervals[i]);
        auto sa_interval_end = (uint64_t) std::stoi(split_sa_intervals[i + 1]);
        const SA_Interval &sa_interval = std::make_pair(sa_interval_start, sa_interval_end);
        sa_intervals.emplace_back(sa_interval);
    }
    return sa_intervals;
}


Site parse_site(const std::string &sites_part_str) {
    Site site;
    for (const auto &pair_i_v: split(sites_part_str, "@")) {
        std::vector<std::string> site_parts = split(pair_i_v, " ");
        if (site_parts.empty())
            continue;

        auto variant_site_marker = (VariantSiteMarker) std::stoi(site_parts[0]);

        std::vector<int> allele;
        for (uint64_t i = 1; i < site_parts.size(); i++) {
            const auto &allele_element = site_parts[i];
            if (!allele_element.empty())
                allele.push_back(std::stoi(allele_element));
        }

        site.emplace_back(VariantSite(variant_site_marker, allele));
    }
    return site;
}


void parse_kmer_index_entry(KmerIndex &kmers, const std::string &line) {
    const std::vector<std::string> &parts = split(line, "|");

    Kmer encoded_kmer = parse_encoded_kmer(parts[0]);
    const auto crosses_marker = parse_crosses_marker_flag(parts[1]);
    if (!crosses_marker)
        kmers.nonvar_kmers.insert(encoded_kmer);

    const auto sa_intervals = parse_sa_intervals(parts[2]);
    if (!sa_intervals.empty())
        kmers.sa_intervals_map[encoded_kmer] = sa_intervals;
    else
        return;

    Sites sites;
    // start at i=4 because of reverse sa intervals
    for (uint64_t i = 4; i < parts.size(); i++) {
        const auto site = parse_site(parts[i]);
        sites.push_back(site);
    }
    kmers.sites_map[encoded_kmer] = sites;
}


KmerIndex lead_kmer_index(const std::string &encoded_kmers_fname) {
    std::ifstream fhandle;
    fhandle.open(encoded_kmers_fname);

    KmerIndex kmers;
    std::string line;
    while (std::getline(fhandle, line)) {
        parse_kmer_index_entry(kmers, line);
    }

    assert(!kmers.sites_map.empty());
    return kmers;
}


KmerIndex get_kmer_index(const std::string &kmer_fname, const PRG_Info &prg_info) {

    const auto encoded_kmers_fname = std::string(kmer_fname) + ".precalc";

    if (!file_exists(encoded_kmers_fname)) {
        std::cout << "Kmer index not found, building..." << std::endl;
        generate_kmer_index(prg_info.allele_mask,
                            kmer_fname,
                            prg_info.max_alphabet_num,
                            prg_info.dna_rank,
                            prg_info.fm_index);
        std::cout << "Finished generating kmer index" << std::endl;
    }

    std::cout << "Loading kmer index from file" << std::endl;
    const auto kmer_index = lead_kmer_index(encoded_kmers_fname);
    return kmer_index;
}
