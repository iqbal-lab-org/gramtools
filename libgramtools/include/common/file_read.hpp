#ifndef GRAMTOOLS_FILE_READ_HPP
#define GRAMTOOLS_FILE_READ_HPP

#include <fstream>

#include <boost/iostreams/filtering_streambuf.hpp>

using namespace boost::iostreams;

static bool is_gzipped(std::string const &fname) {
  if (fname.substr(fname.size() - 2) == "gz") return true;
  return false;
}

filtering_istreambuf &input_fasta(filtering_istreambuf &input_filter,
                                  std::istream &input_stream, bool gzipped);

#endif  // GRAMTOOLS_FILE_READ_HPP
