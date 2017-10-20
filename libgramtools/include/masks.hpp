#include <vector>
#include <string>
#include <fstream>


#ifndef GRAMTOOLS_PARSE_MASKS_H
#define GRAMTOOLS_PARSE_MASKS_H

using SitesMask = std::vector<uint64_t>;
using AlleleMask = std::vector<int>;

class MasksParser {
public:
    uint64_t max_alphabet_num;
    std::vector<uint64_t> sites;

    std::vector<int> allele;
    std::vector<std::vector<double>> allele_coverage;

    MasksParser(){};
    MasksParser(const std::string &sites_fname, const std::string &alleles_fname);

private:
    void parse_sites(std::istream &stream);
    void parse_allele(std::istream &stream);
};

#endif //GRAMTOOLS_PARSE_MASKS_H
