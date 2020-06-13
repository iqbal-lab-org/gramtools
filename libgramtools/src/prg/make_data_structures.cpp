#include "prg/make_data_structures.hpp"
#include <filesystem>
#include "prg/coverage_graph.hpp"

namespace fs = std::filesystem;

using namespace gram;

FM_Index gram::generate_fm_index(BuildParams const &parameters) {
  FM_Index fm_index;

  sdsl::memory_monitor::start();

  sdsl::cache_config config;
  auto construction_tmp_dir =
      fs::path(parameters.gram_dirpath) / fs::path("sdsl_tmp");
  config.dir = fs::absolute(construction_tmp_dir).string();
  fs::create_directories(config.dir);

  // Last param is the number of bytes per integer for reading encoded PRG
  // string. NB: sdsl doc says reads those in big endian, but actually reads in
  // little endian (GH issue #418) So the prg file needs to be in little endian.
  sdsl::construct(fm_index, parameters.encoded_prg_fpath, config,
                  gram::num_bytes_per_integer);
  fs::remove_all(config.dir);
  sdsl::memory_monitor::stop();

  std::ofstream memory_log_fhandle(parameters.sdsl_memory_log_fpath);
  sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(memory_log_fhandle);

  sdsl::store_to_file(fm_index, parameters.fm_index_fpath);
  return fm_index;
}

FM_Index gram::load_fm_index(CommonParameters const &parameters) {
  FM_Index fm_index;
  sdsl::load_from_file(fm_index, parameters.fm_index_fpath);
  return fm_index;
}

coverage_Graph gram::generate_cov_graph(CommonParameters const &parameters,
                                        PRG_String const &prg_string) {
  coverage_Graph c_g{prg_string};

  // Serialise the cov graph
  std::ofstream ofs{parameters.cov_graph_fpath};
  boost::archive::binary_oarchive oa{ofs};
  oa << c_g;

  return c_g;
}

child_map gram::build_child_map(parental_map const &par_map) {
  child_map result;

  Marker child_marker, parental_marker;
  AlleleId parental_haplotype;
  for (auto const &entry : par_map) {
    child_marker = entry.first;
    // Parental locus: pair of (siteID, AlleleId)
    parental_marker = entry.second.first;
    parental_haplotype = entry.second.second;
    assert(parental_haplotype >= FIRST_ALLELE);

    result[parental_marker][parental_haplotype].push_back(child_marker);
  }
  return result;
}

/**************
 * Bit masks **
 **************/

/**
 * Generates a bit vector with bit set if the given DNA `base` is present at
 * each index of the BWT.
 */
void populate_dna_bwt_masks(FM_Index const &fm_index, DNA_BWT_Masks &d_m) {
  for (uint64_t i = 0; i < fm_index.bwt.size(); i++) {
    switch (fm_index.bwt[i]) {
      case 1:
        d_m.mask_a[i] = true;
        break;
      case 2:
        d_m.mask_c[i] = true;
        break;
      case 3:
        d_m.mask_g[i] = true;
        break;
      case 4:
        d_m.mask_t[i] = true;
        break;
    }
  }
}

/**
 * Generates a filename for a BWT mask for nucleotide bases.
 * @see generate_base_bwt_mask()
 */
std::string bwt_mask_fname(const std::string &base_char,
                           CommonParameters const &parameters) {
  auto handling_unit_tests = parameters.gram_dirpath[0] == '@';
  if (handling_unit_tests) {
    return parameters.gram_dirpath + "_" + base_char + "_base_bwt_mask";
  }
  fs::path dir(parameters.gram_dirpath);
  fs::path file(base_char + "_base_bwt_mask");
  fs::path full_path = dir / file;
  return full_path.string();
}

DNA_BWT_Masks gram::generate_bwt_masks(FM_Index const &fm_index,
                                       CommonParameters const &parameters) {
  auto bwt_size = fm_index.bwt.size();
  DNA_BWT_Masks d_m;
  // Set storage for the bit masks
  d_m.mask_a = sdsl::bit_vector(bwt_size, 0);
  d_m.mask_c = sdsl::bit_vector(bwt_size, 0);
  d_m.mask_g = sdsl::bit_vector(bwt_size, 0);
  d_m.mask_t = sdsl::bit_vector(bwt_size, 0);

  populate_dna_bwt_masks(fm_index, d_m);

  auto fpath = bwt_mask_fname("a", parameters);
  sdsl::store_to_file(d_m.mask_a, fpath);

  fpath = bwt_mask_fname("c", parameters);
  sdsl::store_to_file(d_m.mask_c, fpath);

  fpath = bwt_mask_fname("g", parameters);
  sdsl::store_to_file(d_m.mask_g, fpath);

  fpath = bwt_mask_fname("t", parameters);
  sdsl::store_to_file(d_m.mask_t, fpath);

  return d_m;
}

sdsl::bit_vector load_base_bwt_mask(const std::string &base_char,
                                    CommonParameters const &parameters) {
  auto fpath = bwt_mask_fname(base_char, parameters);
  sdsl::bit_vector mask;
  sdsl::load_from_file(mask, fpath);
  return mask;
}

DNA_BWT_Masks gram::load_dna_bwt_masks(const FM_Index &fm_index,
                                       CommonParameters const &parameters) {
  DNA_BWT_Masks dna_bwt_masks;
  dna_bwt_masks.mask_a = load_base_bwt_mask("a", parameters);
  dna_bwt_masks.mask_c = load_base_bwt_mask("c", parameters);
  dna_bwt_masks.mask_g = load_base_bwt_mask("g", parameters);
  dna_bwt_masks.mask_t = load_base_bwt_mask("t", parameters);
  return dna_bwt_masks;
}

sdsl::bit_vector gram::generate_bwt_markers_mask(const FM_Index &fm_index) {
  sdsl::bit_vector bwt_markers_mask(fm_index.bwt.size(), 0);
  for (uint64_t i = 0; i < fm_index.bwt.size(); i++)
    bwt_markers_mask[i] = fm_index.bwt[i] > 4;
  return bwt_markers_mask;
}
