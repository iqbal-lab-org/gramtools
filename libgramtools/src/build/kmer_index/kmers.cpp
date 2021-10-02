#include "build/kmer_index/kmers.hpp"

#include <build/parameters.hpp>

using namespace gram;

std::vector<Sequence> gram::reverse(
    const ordered_vector_set<Sequence> &reverse_kmers) {
  std::vector<Sequence> kmers;
  for (auto reverse_kmer : reverse_kmers) {
    std::reverse(reverse_kmer.begin(), reverse_kmer.end());
    auto &kmer = reverse_kmer;
    kmers.emplace_back(kmer);
  }
  return kmers;
}

/**
 * Given the current pattern, find the next one.
 * The rightmost incrementable (value < 4) position is the one incremented,
 * maximising prefix conservation.
 */
void next_kmer(Sequence &current_kmer, const uint64_t &kmer_size) {
  int64_t max_update_index = kmer_size - 1;

  // TODO: memory leakage here. Replace with: while (max_update_index >= 0 and
  // current_kmer[max_update_index] == 4)
  while (current_kmer[max_update_index] == 4) max_update_index--;

  if (max_update_index < 0) {  // We have reached '4 4 4 4' and so we are done
    current_kmer = {};
    return;
  }
  // Increment the focal position
  current_kmer[max_update_index] = current_kmer[max_update_index] + (uint8_t)1;
  // Reset to 1 all positions to the right of the incremented position
  for (uint64_t i = (uint64_t)max_update_index + 1; i < kmer_size; i++)
    current_kmer[i] = 1;
}

std::vector<Sequence> gram::get_prefix_diffs(
    const std::vector<Sequence> &kmers) {
  std::vector<Sequence> prefix_diffs = {};
  Sequence last_full_kmer = {};

  for (const auto &kmer : kmers) {
    if (last_full_kmer.empty()) {
      last_full_kmer = kmer;
      prefix_diffs.push_back(last_full_kmer);
      continue;
    }

    bool prefix_found_flag = false;
    std::list<int_Base> prefix_diff_list = {};

    for (int64_t i = last_full_kmer.size() - 1; i >= 0; --i) {
      auto &base = kmer[i];  // The current base in the current kmer
      auto &last_full_base =
          last_full_kmer[i];  // The current base in the previous kmer

      if (base != last_full_base)  // Found first difference from the right
        prefix_found_flag = true;

      if (prefix_found_flag) prefix_diff_list.push_front(base);
    }
    last_full_kmer = kmer;  // Update the kmer predecessor
    // Produce a `Pattern` (vector of `Base`s) from the list of differences
    Sequence prefix_diff{std::make_move_iterator(std::begin(prefix_diff_list)),
                         std::make_move_iterator(std::end(prefix_diff_list))};
    prefix_diffs.push_back(prefix_diff);
  }
  return prefix_diffs;
}

ordered_vector_set<Sequence> gram::generate_all_kmers(
    const uint64_t &kmer_size) {
  ordered_vector_set<Sequence> all_kmers = {};
  Sequence current_kmer(kmer_size, 1);  // Start with the pattern '1 1 1 1'

  while (true) {
    all_kmers.insert(current_kmer);
    next_kmer(current_kmer, kmer_size);
    if (current_kmer.empty()) break;
  }
  return all_kmers;
}

std::vector<Sequence> gram::get_all_kmers(const uint64_t &kmer_size) {
  auto ordered_reverse_kmers = generate_all_kmers(kmer_size);
  // Call to reverse: changes for eg '1234' to '4321'. c[j]=c[kmers_size-i-1], i
  // the original position, j the new. Then the kmers are stored as seen in the
  // prg, but in ordered fashion such that they have maximally identical
  // suffixes.
  auto kmers = reverse(ordered_reverse_kmers);
  return kmers;
}

std::vector<Sequence> gram::get_all_kmer_and_compute_prefix_diffs(
    BuildParams const &parameters, const PRG_Info &prg_info) {
  std::cout << "Getting all kmers" << std::endl;
  auto kmers = get_all_kmers(parameters.kmers_size);
  std::cout << "Getting kmer prefix diffs" << std::endl;
  auto prefix_diffs = get_prefix_diffs(kmers);
  return prefix_diffs;
}
