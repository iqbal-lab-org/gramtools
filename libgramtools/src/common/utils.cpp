#include <cstdint>
#include <vector>
#include <string>
#include <iostream>

#include <boost/filesystem.hpp>

#include "sequence_read/seqread.hpp"
#include "common/utils.hpp"

namespace fs = boost::filesystem;
using namespace gram;

std::string gram::mkdir(std::string const& parent_dirpath, std::string const& child_dirpath){
    fs::path dir1(parent_dirpath), dir2(child_dirpath);
    assert(fs::exists(dir1));
    fs::path full_path = dir1 / dir2;
    if (! fs::exists(full_path)) fs::create_directory(full_path);
    return full_path.string();
}


/******************
 * Data typedefs **
 ******************/

child_map gram::build_child_map(parental_map const& par_map){
    child_map result;

    Marker child_marker, parental_marker;
    AlleleId parental_haplotype;
    for (auto const& entry : par_map){
        child_marker = entry.first;
        // Parental locus: pair of (siteID, AlleleId)
        parental_marker = entry.second.first;
        parental_haplotype = entry.second.second;
        // 1-based in par_map and we move to 0-based in child_map
        assert(parental_haplotype >= 1);

        result[parental_marker][parental_haplotype - 1].push_back(child_marker);
    }
    return result;
}


bool gram::is_site_marker(Marker const& variant_marker){
    if (!(variant_marker > 4)) throw std::invalid_argument("The given marker is not a variant marker (>4)");
    return variant_marker % 2 == 1;
}

bool gram::is_allele_marker(Marker const& variant_marker){
    return !is_site_marker(variant_marker);
}

void gram::ensure_is_site_marker(Marker const& site_ID){
    if (!is_site_marker(site_ID)) throw std::invalid_argument("The given marker is not a site ID");
}

std::size_t gram::siteID_to_index(Marker const& site_ID){
    ensure_is_site_marker(site_ID);
    return (site_ID - 5) / 2;
};

Marker gram::index_to_siteID(std::size_t const& idx){
    return idx * 2 + 5;
}

/******************************************
 * Characters to integers and vice-versa **
 ******************************************/
EncodeResult gram::encode_char(const char &c) {
    EncodeResult encode_result = {};

    switch (c) {
        case 'A':
        case 'a':
            encode_result.is_dna = true;
            encode_result.character = 1;
            return encode_result;

        case 'C':
        case 'c':
            encode_result.is_dna = true;
            encode_result.character = 2;
            return encode_result;

        case 'G':
        case 'g':
            encode_result.is_dna = true;
            encode_result.character = 3;
            return encode_result;

        case 'T':
        case 't':
            encode_result.is_dna = true;
            encode_result.character = 4;
            return encode_result;

        default:
            /// The character is non-DNA, so must be a variant marker.
            encode_result.is_dna = false;
            encode_result.character = (uint32_t) c - '0';
            return encode_result;
    }
}

int_Base gram::encode_dna_base(const char &base_str) {
    EncodeResult result = encode_char(base_str);
    if (result.is_dna) return result.character;
    else return 0;
}

std::string  gram::decode_dna_base(const int_Base& base){
    switch(base){
        case 1: return "A";
        case 2: return "C";
        case 3: return "G";
        case 4: return "T";
        default: std::cerr << "Error: the argument " << base << " is not in [1,4]"; exit(1);
    }
}

Sequence gram::encode_dna_bases(const std::string &dna_str) {
    Sequence pattern;
    for (const auto &base_str: dna_str) {
        int_Base encoded_base = encode_dna_base(base_str);
        if (encoded_base == 0)
            return Sequence {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}


Sequence gram::encode_dna_bases(const GenomicRead &read_sequence) {
    const auto sequence_length = read_sequence.seq.size();
    Sequence pattern;
    for (uint64_t i = 0; i < sequence_length; i++) {
        int_Base encoded_base = encode_dna_base(read_sequence.seq[i]);
        if (encoded_base == 0)
            return Sequence {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}


/************************************************
 * Linearised PRGs to integer vectors, and v-v.**
 ************************************************/
std::string gram::ints_to_prg_string(std::vector<Marker> const& int_vec){
    std::string readable_string(int_vec.size(), '0');
    std::unordered_map<int,int> last_allele_indices; // Will record where to close the sites.

    int pos{-1};
    for (auto& s : int_vec){
        pos++;
        if (s > 4){
            if (s%2 == 1) readable_string[pos] = '[';
            else{
                readable_string[pos] = ',';
                if (last_allele_indices.find(s) != last_allele_indices.end()){
                    last_allele_indices.erase(s);
                }
                last_allele_indices.insert({s, pos});
            }
            continue;
        }
        // Implicit else: s <= 4
        char base = decode_dna_base(s).c_str()[0];
        readable_string[pos] = base;
    }

    // Close the sites
    for (auto s : last_allele_indices){
        auto pos = s.second;
        readable_string[pos] = ']';
    }
    return readable_string;
}

marker_vec gram::prg_string_to_ints(std::string const& string_prg) {
    int_Base base;
    std::stack<int> marker_stack;
    int max_var_marker{3};
    int char_count{0};

    marker_vec encoded_prg(string_prg.size());
    for (std::size_t i = 0; i < string_prg.size(); ++i){
        const auto &c = string_prg[i];

        switch(c) {
            case '[' : {
                max_var_marker += 2;
                marker_stack.push(max_var_marker);
                encoded_prg[char_count++] = max_var_marker;
                break;
            }

            case ']' : {
                assert(!marker_stack.empty());
                encoded_prg[char_count++] = marker_stack.top() + 1;
                marker_stack.pop();
                break;
            }

            case ',' : {
                assert(!marker_stack.empty());
                encoded_prg[char_count++] = marker_stack.top() + 1;
                break;
            }

            default : {
                    base = encode_dna_base(c);
                    if (base == 0){
                        std::cerr << "Error: the argument " << c << " is not a nucleotide char";
                        exit(1);
                    }
                    encoded_prg[char_count++] = encode_dna_base(c);
                    break;
                }
            }
        }

    // BOOST_LOG_TRIVIAL(info) << "Number of sites produced: " << (max_var_marker -3 ) / 2;
    return encoded_prg;
}
