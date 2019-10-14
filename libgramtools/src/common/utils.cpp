#include <cstdint>
#include <vector>
#include <string>
#include <iostream>

#include <boost/filesystem.hpp>

#include "sequence_read/seqread.hpp"
#include "common/utils.hpp"


namespace fs = boost::filesystem;
using namespace gram;


bool gram::is_site_marker(Marker const& variant_marker){
    if (!variant_marker > 4) throw std::invalid_argument("The given marker is not a variant marker");
    return variant_marker % 2 == 1;
}

bool gram::is_allele_marker(Marker const& variant_marker){
    return !is_site_marker(variant_marker);
}

std::string gram::full_path(const std::string &gram_dirpath,
                      const std::string &file_name) {
    fs::path dir(gram_dirpath);
    fs::path file(file_name);
    fs::path full_path = dir / file;
    return full_path.string();
}


int_Base gram::encode_dna_base(const char &base_str) {
    switch (base_str) {
        case 'A':
        case 'a': return 1;

        case 'C':
        case 'c': return 2;

        case 'G':
        case 'g': return 3;

        case 'T':
        case 't': return 4;

        default: return 0;
    }
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

Pattern gram::encode_dna_bases(const std::string &dna_str) {
    Pattern pattern;
    for (const auto &base_str: dna_str) {
        int_Base encoded_base = encode_dna_base(base_str);
        if (encoded_base == 0)
            return Pattern {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}


/**
 * Produce integer-encoded Watson-Crick base complement.
 */
int_Base complement_encoded_base(const int_Base &encoded_base) {
    switch (encoded_base) {
        case 1:
            return 4;
        case 2:
            return 3;
        case 3:
            return 2;
        case 4:
            return 1;
        default:
            return 0;
    }
}


Pattern gram::reverse_complement_read(const Pattern &read) {
    Pattern reverse_read;
    reverse_read.reserve(read.size());

    for (auto it = read.rbegin(); it != read.rend(); ++it) {
        const auto &base = *it;
        auto compliment_base = complement_encoded_base(base);
        reverse_read.push_back(compliment_base);
    }
    return reverse_read;
}


Pattern gram::encode_dna_bases(const GenomicRead &read_sequence) {
    const auto sequence_length = strlen(read_sequence.seq);
    Pattern pattern;
    for (uint64_t i = 0; i < sequence_length; i++) {
        int_Base encoded_base = encode_dna_base(read_sequence.seq[i]);
        if (encoded_base == 0)
            return Pattern {};
        pattern.emplace_back(encoded_base);
    }
    return pattern;
}


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

std::vector<Marker> gram::prg_string_to_ints(std::string const& string_prg) {
    int_Base base;
    std::stack<int> marker_stack;
    int max_var_marker{3};
    int char_count{0};

    std::vector<Marker> encoded_prg(string_prg.size());
    for (int i = 0; i<string_prg.size(); ++i){
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
