#include "common/file_read.hpp"

#include <boost/iostreams/categories.hpp> // input_filter_tag
#include <boost/iostreams/operations.hpp> // get, WOULD_BLOCK
#include <boost/iostreams/filter/gzip.hpp> // decompressor

struct to_upper_filter {
    typedef char              char_type;
    typedef input_filter_tag  category;

    template<typename Source>
    int get(Source& src)
    {
        int c = boost::iostreams::get(src);
        if (c == EOF) return c;
        return std::toupper((unsigned char) c);
    }
};

filtering_istreambuf &input_fasta(filtering_istreambuf &input_filter, std::istream &input_stream, bool gzipped) {
    input_filter.push(to_upper_filter());
    if (gzipped){
        input_filter.push(gzip_decompressor());
    }
    input_filter.push(input_stream);
    return input_filter;
}

