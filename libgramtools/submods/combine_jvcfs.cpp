/**
 * @file Combine JSON genotyped files into one
 */
#include <filesystem>
#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

#include "genotype/infer/output_specs/json_prg_spec.hpp"

namespace fs = std::filesystem;
using JSON = nlohmann::json;
using namespace gram::json;

void usage(const char* argv[]) {
  std::cout << "Usage: " << argv[0] << " fofn fout" << std::endl;
  std::cout << "\t fofn: file of file names of the JSON files to combine"
            << std::endl;
  std::cout << "\t fout: name of output combined JSON file" << std::endl;
  exit(1);
}

int main(int argc, const char* argv[]) {
  if (argc != 3) usage(argv);
  fs::path fofn(argv[1]);
  if (!fs::exists(fofn)) {
    std::cout << fofn << " not found.";
    usage(argv);
  } else if (fs::is_empty(fofn)) {
    std::cout << fofn << " is empty.";
    usage(argv);
  }

  std::ofstream fout(argv[2]);
  if (!fout.good()) {
    std::cout << "Error: could not open " << argv[2] << std::endl;
    usage(argv);
  }

  std::ifstream fin(fofn);
  std::string next_file;
  JSON next_json;
  Json_Prg combined_json;
  bool first{true};
  while (std::getline(fin, next_file)) {
    std::ifstream next_istream(next_file);
    if (!next_istream.good()) {
      std::cout << "Error: Could not open JSON file " << next_file << std::endl;
      exit(1);
    }
    next_json.clear();
    next_istream >> next_json;
    if (first) {
      combined_json.set_prg(next_json);
      first = false;
    }

    else {
      Json_Prg next_prg(next_json);
      combined_json.combine_with(next_prg);
    }
  }
  fout << combined_json.get_prg();
}
