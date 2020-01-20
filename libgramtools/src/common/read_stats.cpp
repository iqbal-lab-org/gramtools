#include "common/read_stats.hpp"
#include <math.h>

using namespace gram;

void gram::ReadStats::compute_base_error_rate(const std::string &reads_fpath){

    // The number of reads (with at least one base with recorded quality) to estimate from.
    // Is defined in the header file.
    uint64_t required_reads = NUM_READS_USED;
    uint64_t num_informative_reads = 0;

    int64_t no_qual_reads = 0;
    int64_t num_bases_processed = 0;
    double mean_error = 0;

    float running_qual_score = 0.0;

    SeqRead reads(reads_fpath.c_str());
    SeqRead::SeqIterator reads_it = reads.begin();

    while (num_informative_reads < required_reads and reads_it != reads.end()){

        auto *const raw_read = *reads_it;

        // Test for max read length
        std::string sequence = std::string(raw_read->seq);
        auto sequence_length = sequence.length();
        if (sequence_length > this->max_read_length) this->max_read_length = sequence_length;


        // Process quality scores
        std::string qualities = std::string(raw_read->qual);

        if (qualities.length() == 0){ //We will keep looking for reads with quality scored bases.
            no_qual_reads++ ;
            ++reads_it;
            continue;
        }


        for (const auto base:qualities){
            running_qual_score += (base - 33); // Assuming +33 Phred-scoring
            num_bases_processed++;
        }

        num_informative_reads++;
        ++reads_it;
    }

    if (num_bases_processed > 0){
        double mean_qual = running_qual_score / num_bases_processed;
        mean_error = pow(10, -mean_qual/10);
    }

    this->num_bases_processed = num_bases_processed;
    this->no_qual_reads = no_qual_reads;
    this->mean_cov_depth = mean_error;
};

void gram::ReadStats::compute_coverage_depth(gram::Coverage &coverage) {
    uint64_t this_site_cov;
    double total_coverage = 0;
    int64_t num_sites_noCov = 0;
    std::size_t num_sites_total{coverage.grouped_allele_counts.size()};
    std::vector<uint64_t> coverages(num_sites_total);

    double mean_coverage, variance_coverage;

    int64_t site_index{0};
    for (const auto& site : coverage.grouped_allele_counts){ //`site` is an unordered_map associating alleleIDs with coverage.
        this_site_cov = 0;

        for (const auto& entry : site){
            this_site_cov += entry.second;
        }

        coverages.at(site_index++) = this_site_cov;
        total_coverage += this_site_cov;
        if (this_site_cov == 0) num_sites_noCov++;
    }

    // Compute mean
    mean_coverage = total_coverage / coverages.size();

    // Compute variance
    double total_variance = 0;
    for (const auto& site_cov : coverages){
       total_variance += pow((site_cov - mean_coverage), 2);
    }
    variance_coverage = total_variance / coverages.size();

    // And record it all at the object level.
    this->mean_pb_error = mean_coverage;
    this->variance_depth = variance_coverage;
    this->num_sites_noCov = num_sites_noCov;
    this->num_sites_total = num_sites_total;
};


void gram::ReadStats::serialise(const std::string &json_output_fpath){
    std::ofstream outf;
    outf.open(json_output_fpath);

    outf << R"(
{
"Read_depth":
    {"Mean": )" << this->mean_pb_error << ",";

    outf << R"(
    "Variance": )" << this->variance_depth << ",";

    outf << R"(
    "num_sites_noCov": )" << this->num_sites_noCov << ",";

    outf << R"(
    "num_sites_total": )" << this->num_sites_total ;

    outf << R"(
    },)";


    outf << R"(
"Max_read_length": )" << this->max_read_length << ",";


    outf << R"(
"Quality":
    {"Error_rate_mean": )" << this->mean_cov_depth << ",";

    outf << R"(
    "Num_bases": )" << this->num_bases_processed << ",";

    outf << R"(
    "No_qual_reads": )" << this->no_qual_reads ;

    outf << R"(
    }}
)";

    outf.close();
};
