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


void dump_kmer_index(std::ofstream &kmer_index_file, const KmerIndex &kmer_index) {
    for (const auto &kmer_sa_intervals: kmer_index.sa_intervals_map) {
        const auto &kmer = kmer_sa_intervals.first;
        const SA_Intervals &sa_intervals = kmer_index.sa_intervals_map.at(kmer);
        const std::string kmer_entry =
                dump_kmer_index_entry(kmer,
                                      sa_intervals,
                                      kmer_index.nonvar_kmers,
                                      kmer_index.sites_map);
        kmer_index_file << kmer_entry << std::endl;
    }
}


void print_sa_intervals(const SA_Intervals &sa_intervals) {
    std::cout << "sa_intervals size: " << sa_intervals.size() << std::endl;
    for (const auto &sa_interval: sa_intervals)
        std::cout << "(" << sa_interval.first << ", " << sa_interval.second << ")" << "   ";
    if (sa_intervals.size() > 0)
        std::cout << std::endl;
}


void print_sites(const Sites &sites) {
    std::cout << "sites length: " << sites.size() << std::endl;
    for (const auto &site: sites) {
        std::cout << "site length: " << site.size() << std::endl;
        for (const auto &variant_site: site) {
            std::cout << "variant site marker: " << variant_site.first << std::endl;
            for (const auto &allele: variant_site.second)
                std::cout << allele << ", ";
            if (variant_site.second.size() > 0)
                std::cout << std::endl;
        }
    }
}


void print_cache(const KmerIndexCache &cache) {
    for (auto &elem: cache) {
        std::cout << "***start elem" << std::endl;
        std::cout << (int) elem.base << std::endl;

        print_sa_intervals(elem.sa_intervals);
        print_sites(elem.sites);
        std::cout << "***end elem" << std::endl;
        std::cout << std::endl;
    }
}


CacheElement get_next_cache_element(SA_Intervals &sa_intervals,
                                    Sites &sites,
                                    const uint8_t base,
                                    const PRG_Info &prg_info) {
    CacheElement new_cache_element;
    new_cache_element.sa_intervals = sa_intervals;
    new_cache_element.sites = sites;

    bool delete_first_interval = false;
    const bool kmer_index_done = false;
    // const bool last_base = it + 1 == kmer_suffix_diff.rend();
    const bool last_base = false;

    reduce_search_scope(base,
                        new_cache_element.sa_intervals,
                        new_cache_element.sites,
                        delete_first_interval,
                        kmer_index_done,
                        last_base,
                        prg_info);

    new_cache_element.base = base;
    return new_cache_element;
}


KmerIndexCache initial_kmer_index_cache(const Kmer &full_kmer,
                                        const PRG_Info &prg_info) {
    KmerIndexCache cache;

    // reverse iteration normally provided by bidir_search_bwd, but we skip that func
    for (auto it = full_kmer.rbegin(); it != full_kmer.rend(); ++it) {
        const auto &base = *it;

        if (cache.empty()) {
            SA_Intervals sa_intervals = {{0, prg_info.fm_index.size()}};
            Sites sites = {Site()};
            auto new_cache_element = get_next_cache_element(sa_intervals, sites,
                                                            base, prg_info);
            cache.emplace_back(new_cache_element);
            continue;
        }

        auto &last_element = cache.back();
        const auto &new_cache_element = get_next_cache_element(last_element.sa_intervals,
                                                               last_element.sites,
                                                               base, prg_info);
        cache.emplace_back(new_cache_element);
    }
    return cache;
}


void update_kmer_index_cache(KmerIndexCache &cache,
                             const Kmer &kmer_suffix_diff,
                             const int kmer_size,
                             const PRG_Info &prg_info) {

    if (kmer_suffix_diff.size() == kmer_size) {
        auto &full_kmer = kmer_suffix_diff;
        cache = initial_kmer_index_cache(full_kmer, prg_info);
        return;
    }

    // truncate cache
    const auto new_cache_size = kmer_size - kmer_suffix_diff.size();
    cache.resize(new_cache_size);

    /*
    auto start_erase = cache.begin();
    std::advance(start_erase, new_cache_size);
    cache.erase(start_erase, ++cache.end());
     */

    // reverse iteration normally provided by bidir_search_bwd, but we skip that func
    for (auto it = kmer_suffix_diff.rbegin(); it != kmer_suffix_diff.rend(); ++it) {
        const auto &base = *it;

        auto &last_element = cache.back();
        const auto &new_cache_element = get_next_cache_element(last_element.sa_intervals,
                                                               last_element.sites,
                                                               base, prg_info);
        cache.emplace_back(new_cache_element);
    }
}


void update_full_kmer(Kmer &full_kmer, const Kmer &kmer_suffix_diff, const int kmer_size) {
    if (kmer_suffix_diff.size() == kmer_size) {
        full_kmer = kmer_suffix_diff;
        return;
    }

    auto start_idx = 0;
    for (const auto &base: kmer_suffix_diff)
        full_kmer[start_idx++] = base;
}


KmerIndex index_kmers(const Kmers &kmer_suffix_diffs,
                      const int kmer_size,
                      const PRG_Info &prg_info) {
    KmerIndex kmer_index;
    Kmer full_kmer;
    KmerIndexCache cache;

    for (auto &kmer_suffix_diff: kmer_suffix_diffs) {
        update_full_kmer(full_kmer,
                         kmer_suffix_diff,
                         kmer_size);

        update_kmer_index_cache(cache,
                                kmer_suffix_diff,
                                kmer_size,
                                prg_info);

        const auto &last_cache_element = cache.back();
        kmer_index.sa_intervals_map[full_kmer] = last_cache_element.sa_intervals;
        kmer_index.sites_map[full_kmer] = last_cache_element.sites;

        if (last_cache_element.sa_intervals.empty()) {
            kmer_index.sa_intervals_map.erase(full_kmer);
        }

        const auto &first_sties_element = last_cache_element.sites.front();
        if (first_sties_element.empty())
            kmer_index.nonvar_kmers.insert(full_kmer);
    }
    return kmer_index;
}


void generate_kmer_index(const std::string &kmer_fname, const int kmer_size, const PRG_Info &prg_info) {
    std::ifstream kmer_fhandle;
    kmer_fhandle.open(kmer_fname);

    Kmers kmer_suffix_diffs;
    std::string line;
    while (std::getline(kmer_fhandle, line)) {
        const Kmer &kmer_suffix_diff = encode_dna_bases(line);
        kmer_suffix_diffs.emplace_back(kmer_suffix_diff);
    }

    KmerIndex kmer_index = index_kmers(kmer_suffix_diffs, kmer_size, prg_info);
    std::ofstream kmer_index_file;
    kmer_index_file.open(std::string(kmer_fname) + ".precalc");
    dump_kmer_index(kmer_index_file, kmer_index);
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


KmerIndex load_kmer_index(const std::string &encoded_kmers_fname) {
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


KmerIndex get_kmer_index(const std::string &kmer_fname, const int kmer_size, const PRG_Info &prg_info) {

    const auto encoded_kmers_fname = std::string(kmer_fname) + ".precalc";

    if (!file_exists(encoded_kmers_fname)) {
        std::cout << "Kmer index not found, building..." << std::endl;
        generate_kmer_index(kmer_fname, kmer_size, prg_info);
        std::cout << "Finished generating kmer index" << std::endl;
    }

    std::cout << "Loading kmer index from file" << std::endl;
    const auto kmer_index = load_kmer_index(encoded_kmers_fname);
    return kmer_index;
}
