#include "genotype/read_stats.hpp"

#include <math.h>

#include "genotype/infer/types.hpp"
#include "prg/coverage_graph.hpp"

using namespace gram;

void gram::ReadStats::compute_base_error_rate(const std::string& reads_fpath) {
  SeqRead reads(reads_fpath.c_str());
  SeqRead::SeqIterator reads_it = reads.begin();
  process_read_perbase_error_rates(reads_it);
}

void gram::ReadStats::compute_base_error_rate(GenomicRead_vector const& reads) {
  GenomicReadIterator reads_it(reads);
  process_read_perbase_error_rates(reads_it);
}

void gram::ReadStats::process_read_perbase_error_rates(
    AbstractGenomicReadIterator& reads_it) {
  // The number of reads (with at least one base with recorded quality) to
  // estimate from. Is defined in the header file.
  uint64_t required_reads = NUM_READS_USED;
  uint64_t num_informative_reads = 0;

  int64_t no_qual_reads = 0;
  int64_t num_bases_processed = 0;
  double mean_error = 0;

  float running_qual_score = 0.0;

  while (num_informative_reads < required_reads and reads_it.has_more_reads()) {
    auto* const raw_read = *reads_it;

    // Test for max read length
    std::string sequence = std::string(raw_read->seq);
    auto sequence_length = sequence.length();
    if (sequence_length > this->max_read_length)
      this->max_read_length = sequence_length;

    // Process quality scores
    std::string qualities = std::string(raw_read->qual);

    if (qualities.length() ==
        0) {  // We will keep looking for reads with quality scored bases.
      no_qual_reads++;
      ++reads_it;
      continue;
    }

    for (const auto base : qualities) {
      running_qual_score += (base - 33);  // Assuming +33 Phred-scoring
      num_bases_processed++;
    }

    num_informative_reads++;
    ++reads_it;
  }

  if (num_bases_processed > 0) {
    double mean_qual = running_qual_score / num_bases_processed;
    mean_error = pow(10, -mean_qual / 10);
  }

  this->num_bases_processed = num_bases_processed;
  this->no_qual_reads = no_qual_reads;
  this->mean_pb_error = mean_error;
}

ReadStats::haplogroup_cov ReadStats::get_max_cov_haplogroup(
    GroupedAlleleCounts const& gped_cov) {
  std::map<AlleleId, CovCount> counts;
  for (auto const& entry : gped_cov) {
    for (auto const& allele_id : entry.first) {
      if (counts.find(allele_id) == counts.end())
        counts.insert(haplogroup_cov{allele_id, 0});
      counts.at(allele_id) += entry.second;
    }
  }

  auto max_elem = std::max_element(
      counts.begin(), counts.end(),
      [](haplogroup_cov const& elem1, haplogroup_cov const& elem2) {
        return elem1.second < elem2.second;
      });
  if (max_elem == counts.end())
    return haplogroup_cov{0, 0};
  else
    return *max_elem;
}

ReadStats::allele_and_cov ReadStats::extract_max_coverage_allele(
    SitesGroupedAlleleCounts const& gped_covs, covG_ptr start_node,
    covG_ptr end_node) {
  Allele result;
  covG_ptr cur_Node = start_node;
  auto site_index = siteID_to_index(cur_Node->get_site_ID());
  auto max_elem = get_max_cov_haplogroup(gped_covs.at(site_index));
  auto const allele_cov = max_elem.second;

  while (cur_Node != end_node) {
    if (cur_Node->is_bubble_start()) {
      cur_Node = cur_Node->get_edges().at(max_elem.first);
      site_index = siteID_to_index(cur_Node->get_site_ID());
      max_elem = get_max_cov_haplogroup(gped_covs.at(site_index));
      continue;
    }
    if (cur_Node->has_sequence()) {
      result =
          result + Allele{cur_Node->get_sequence(), cur_Node->get_coverage()};
    }
    cur_Node = cur_Node->get_edges().at(0);
  }
  return allele_and_cov{result, allele_cov};
}

void gram::AbstractReadStats::compute_coverage_depth(
    Coverage const& coverage, coverage_Graph const& cov_graph) {
  ReadStats::allele_and_cov site_extraction;
  double site_perbase_coverage = 0, total_coverage = 0;
  int64_t num_sites_noCov = 0;
  std::vector<double> coverages;

  double mean_coverage, variance_coverage;

  for (const auto& node_pair : cov_graph.bubble_map) {
    auto site_ID = node_pair.first->get_site_ID();

    // If the site is nested within another, we do not process its coverage
    if (cov_graph.par_map.find(site_ID) != cov_graph.par_map.end()) continue;
    site_extraction = extract_max_coverage_allele(
        coverage.grouped_allele_counts, node_pair.first, node_pair.second);

    if (site_extraction.first.pbCov.size() > 0)
      site_perbase_coverage = site_extraction.first.get_average_cov();
    else  // Case: direct deletion allele
      site_perbase_coverage = (double)site_extraction.second;
    total_coverage += site_perbase_coverage;
    coverages.push_back(site_perbase_coverage);
    if (site_extraction.second == 0) num_sites_noCov++;
  }

  // Compute mean
  mean_coverage = total_coverage / coverages.size();

  // Compute variance
  double total_variance = 0;
  for (const auto& site_cov : coverages) {
    total_variance += pow((site_cov - mean_coverage), 2);
  }
  variance_coverage = total_variance / coverages.size();

  // And record it all at the object level.
  this->mean_cov_depth = mean_coverage;
  this->variance_cov_depth = variance_coverage;
  this->num_sites_noCov = num_sites_noCov;
  this->num_sites_total = coverages.size();
}

void gram::ReadStats::serialise(const std::string& json_output_fpath) {
  std::ofstream outf;
  outf.open(json_output_fpath);

  outf << R"(
{
"Read_depth":
    {"Mean": )"
       << this->mean_cov_depth << ",";

  outf << R"(
    "Variance": )"
       << this->variance_cov_depth << ",";

  outf << R"(
    "num_sites_noCov": )"
       << this->num_sites_noCov << ",";

  outf << R"(
    "num_sites_total": )"
       << this->num_sites_total;

  outf << R"(
    },)";

  outf << R"(
"Max_read_length": )"
       << this->max_read_length << ",";

  outf << R"(
"Quality":
    {"Error_rate_mean": )"
       << this->mean_pb_error << ",";

  outf << R"(
    "Num_bases": )"
       << this->num_bases_processed << ",";

  outf << R"(
    "No_qual_reads": )"
       << this->no_qual_reads;

  outf << R"(
    }}
)";

  outf.close();
}
