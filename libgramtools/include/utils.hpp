#include <vector>
#include <list>
#include <cstdint>
#include <string>

#include "sequence_read/seqread.hpp"


#ifndef GRAMTOOLS_UTILS_HPP
#define GRAMTOOLS_UTILS_HPP

using Base = uint8_t;
using Pattern = std::vector<Base>;
using Patterns = std::vector<Pattern>;

using Marker = uint64_t;
using AlleleId = uint64_t;

using VariantSite = std::pair<Marker, AlleleId>;
using VariantSitePath = std::list<VariantSite>;
using VariantSitePaths = std::list<VariantSitePath>;

using SA_Index = uint64_t;
using SA_Interval = std::pair<SA_Index, SA_Index>;

Pattern encode_dna_bases(const std::string &dna_str);
Pattern encode_dna_bases(const GenomicRead &read_sequence);

#endif //GRAMTOOLS_UTILS_HPP
