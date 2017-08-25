#include "common.hpp"


void cleanup_files() {
    std::remove("csa_file");
    std::remove("csa_file_rev");
    std::remove("int_alphabet_file");
    std::remove("int_alphabet_file_rev");
    std::remove("memory_log_file");
}