#include "prg/make_data_structures.hpp"


using namespace gram;


FM_Index gram::load_fm_index(const Parameters &parameters) {
    FM_Index fm_index;
    sdsl::load_from_file(fm_index, parameters.fm_index_fpath);
    return fm_index;
}

FM_Index gram::generate_fm_index(const Parameters &parameters) {
    FM_Index fm_index;

    sdsl::memory_monitor::start();
    // Last param is the number of bytes per integer for reading encoded PRG string.
    // NB: sdsl doc says reads those in big endian, but actually reads in little endian (GH issue #418)
    // So the prg file needs to be in little endian.
    sdsl::construct(fm_index, parameters.encoded_prg_fpath, gram::num_bytes_per_integer);
    sdsl::memory_monitor::stop();

    std::ofstream memory_log_fhandle(parameters.sdsl_memory_log_fpath);
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(memory_log_fhandle);

    sdsl::store_to_file(fm_index, parameters.fm_index_fpath);
    return fm_index;
}

coverage_Graph gram::generate_cov_graph(const Parameters &parameters, PRG_String const &prg_string){
    coverage_Graph c_g{prg_string};

    // Serialise the cov graph
    std::ofstream ofs{parameters.cov_graph_fpath};
    boost::archive::binary_oarchive oa{ofs};
    oa << c_g;

    return c_g;
}

