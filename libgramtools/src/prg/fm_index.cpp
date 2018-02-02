#include <fstream>
#include <cinttypes>
#include <sdsl/suffix_arrays.hpp>

#include "common/parameters.hpp"
#include "prg/fm_index.hpp"


FM_Index load_fm_index(const Parameters &parameters) {
    FM_Index fm_index;
    sdsl::load_from_file(fm_index, parameters.fm_index_fpath);
    return fm_index;
}


FM_Index generate_fm_index(const Parameters &parameters) {
    FM_Index fm_index;

    sdsl::memory_monitor::start();
    sdsl::construct(fm_index, parameters.encoded_prg_fpath, 0);
    sdsl::memory_monitor::stop();

    std::ofstream memory_log_fhandle(parameters.sdsl_memory_log_fpath);
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(memory_log_fhandle);

    sdsl::store_to_file(fm_index, parameters.fm_index_fpath);
    return fm_index;
}
