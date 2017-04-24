#include <stdint.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "parse_masks.h"


MasksParser::MasksParser(const std::string &sites_fname, const std::string &alleles_fname) {
    parse_sites(sites_fname);
    parse_allele(alleles_fname);
}


void MasksParser::parse_sites(const std::string &sites_fname) {
    // TODO: Name d clearly
    uint64_t d;
    uint64_t count_sites = 0;
    std::ifstream fhandle(sites_fname);

    while (fhandle >> d) {
        if (d > count_sites)
            count_sites = d;
        MasksParser::sites.push_back(d);
    }
    fhandle.close();

    // cout<<"Covgs dim"<<covgs.size()<<" "<<covgs.front().size()<<" "<<covgs.back().size()<<endl;
    // no_sites is last odd number in mask_sites, but alphabet size
    // is the even number corresponding to it
    MasksParser::maxx = count_sites + 1;
}


void MasksParser::parse_allele(const std::string &alleles_fname) {
    std::ifstream fhandle(alleles_fname);
    int count_allele = 0;
    int a, i = 0;
    // TODO: Assert at least 2 alleles at each site
    while (fhandle >> a) {
        if (a > count_allele)
            count_allele = a;
        if (a < count_allele && a != 0) {
            MasksParser::allele_coverage.push_back(std::vector<int>(count_allele, 0));
            count_allele = a;
        }
        i++;
        MasksParser::allele.push_back(a);
    }
    fhandle.close();

    if (count_allele > 0)
        MasksParser::allele_coverage.push_back(std::vector<int>(count_allele, 0));
}


/*
uint64_t parse_masks(std::vector<uint64_t> &mask_s, std::vector<int> &mask_a,
                     std::string sites_fname, std::string alleles_fname,
                     std::vector<std::vector<int>> &covgs) {

    uint64_t d, no_sites;
    std::ifstream h1(sites_fname);

    no_sites = 0;
    while (h1 >> d) {
        if (d > no_sites)
            no_sites = d;
        mask_s.push_back(d);
    }
    h1.close();

    std::ifstream h2(alleles_fname);
    int no_alleles = 0;
    int a, i = 0;
    // TODO: Assert at least 2 alleles at each site
    while (h2 >> a) {
        if (a > no_alleles)
            no_alleles = a;
        if (a < no_alleles && a != 0) {
            covgs.push_back(std::vector<int>(no_alleles, 0));
            no_alleles = a;
        }
        i++;
        mask_a.push_back(a);
    }
    h2.close();

    if (no_alleles > 0)
        covgs.push_back(std::vector<int>(no_alleles, 0));

    std::cout << std::endl << mask_s.size() << " " << mask_a.size() << std::endl;

    // cout<<"Covgs dim"<<covgs.size()<<" "<<covgs.front().size()<<" "<<covgs.back().size()<<endl;
    // no_sites is last odd number in mask_sites, but alphabet size is the even number corresponding to it
    return (no_sites + 1);
}
*/
