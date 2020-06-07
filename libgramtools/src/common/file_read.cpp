#include "common/file_read.hpp"

#include <boost/iostreams/categories.hpp> // input_filter_tag
#include <boost/iostreams/operations.hpp> // get, WOULD_BLOCK
#include <boost/iostreams/filter/gzip.hpp> // decompressor

struct fasta_to_upper_filter {
    typedef char              char_type;
    typedef input_filter_tag  category;
    bool newline{true}, header{false};

    template<typename Source>
    int get(Source& src)
    {
        int c = boost::iostreams::get(src);
        if (newline){
            (unsigned char)c == '>' ? header = true : header = false;
            newline = false;
        }
        if ((unsigned char)c == '\n') {
            newline = true;
            header = false;
        }
        if (c == EOF || header) return c;
        return std::toupper((unsigned char) c);
    }
};

filtering_istreambuf &input_fasta(filtering_istreambuf &input_filter, std::istream &input_stream, bool gzipped) {
    input_filter.push(fasta_to_upper_filter());
    if (gzipped){
        input_filter.push(gzip_decompressor());
    }
    input_filter.push(input_stream);
    return input_filter;
}

