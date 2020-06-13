#include "common/parameters.hpp"

using namespace gram;

std::string gram::mkdir(std::string const& parent_dirpath,
                        std::string const& child_dirpath) {
  fs::path dir1(parent_dirpath), dir2(child_dirpath);
  dir1 = fs::absolute(dir1);
  assert(fs::exists(dir1));
  fs::path full_path = dir1 / dir2;
  if (!fs::exists(full_path)) fs::create_directory(full_path);
  return full_path.string();
}

std::string gram::full_path(const std::string& base_dirpath,
                            const std::string& file_name) {
  auto dir = fs::absolute(fs::path(base_dirpath));
  fs::path file(file_name);
  fs::path full_path = dir / file;
  return full_path.string();
}

void gram::fill_common_parameters(CommonParameters& parameters,
                                  std::string const gram_dirpath) {
  parameters.gram_dirpath = gram_dirpath;
  parameters.built_vcf =
      full_path(gram_dirpath, "build.vcf");  // May not exist, if --prg used
  parameters.encoded_prg_fpath = full_path(gram_dirpath, "prg");
  parameters.prg_coords_fpath = full_path(gram_dirpath, "prg_coords.tsv");
  parameters.fm_index_fpath = full_path(gram_dirpath, "fm_index");
  parameters.cov_graph_fpath = full_path(gram_dirpath, "cov_graph");
  parameters.sites_mask_fpath = full_path(gram_dirpath, "variant_site_mask");
  parameters.allele_mask_fpath = full_path(gram_dirpath, "allele_mask");

  parameters.kmer_index_fpath = full_path(gram_dirpath, "kmer_index");
  parameters.kmers_fpath = full_path(gram_dirpath, "kmers");
  parameters.kmers_stats_fpath = full_path(gram_dirpath, "kmers_stats");
  parameters.sa_intervals_fpath = full_path(gram_dirpath, "sa_intervals");
  parameters.paths_fpath = full_path(gram_dirpath, "paths");
}
