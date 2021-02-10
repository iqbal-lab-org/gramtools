#include "genotype/quasimap/search/BWT_search.hpp"

#include <sdsl/suffix_arrays.hpp>

using namespace gram;

uint64_t gram::dna_bwt_rank(const uint64_t &upper_index, const Marker &dna_base,
                            const PRG_Info &prg_info) {
  switch (dna_base) {
    case 1:
      return prg_info.rank_bwt_a(upper_index);
    case 2:
      return prg_info.rank_bwt_c(upper_index);
    case 3:
      return prg_info.rank_bwt_g(upper_index);
    case 4:
      return prg_info.rank_bwt_t(upper_index);
    default:
      return 0;
  }
}

/**
 * Backward search followed by check whether the extended searched pattern maps
 * somewhere in the prg.
 */
std::optional<SearchState> search_fm_index_base_backwards(
    const int_Base &pattern_char, const uint64_t char_first_sa_index,
    const SearchState &search_state, const PRG_Info &prg_info) {
  auto next_sa_interval = base_next_sa_interval(
      pattern_char, char_first_sa_index, search_state.sa_interval, prg_info);
  //  An 'invalid' SA interval (i,j) is defined by i-1=j, which occurs when the
  //  read no longer maps anywhere in the prg.
  auto valid_sa_interval =
      next_sa_interval.first - 1 != next_sa_interval.second;
  if (not valid_sa_interval) return {};

  auto new_search_state = search_state;
  new_search_state.sa_interval.first = next_sa_interval.first;
  new_search_state.sa_interval.second = next_sa_interval.second;
  return new_search_state;
}

SA_Interval gram::base_next_sa_interval(
    const Marker &next_char, const SA_Index &next_char_first_sa_index,
    const SA_Interval &current_sa_interval, const PRG_Info &prg_info) {
  const auto &current_sa_start = current_sa_interval.first;
  const auto &current_sa_end = current_sa_interval.second;

  SA_Index sa_start_offset;
  if (current_sa_start <= 0)
    sa_start_offset = 0;
  else {
    //  TODO: Consider deleting this if-clause, next_char should never be > 4,
    //  it probably never runs
    if (next_char > 4)
      sa_start_offset = prg_info.fm_index.bwt.rank(current_sa_start, next_char);
    else {
      sa_start_offset = dna_bwt_rank(current_sa_start, next_char, prg_info);
    }
  }

  SA_Index sa_end_offset;
  //  TODO: Consider deleting this if-clause, next_char should never be > 4, it
  //  probably never runs
  if (next_char > 4)
    sa_end_offset = prg_info.fm_index.bwt.rank(current_sa_end + 1, next_char);
  else {
    sa_end_offset = dna_bwt_rank(current_sa_end + 1, next_char, prg_info);
  }

  auto new_start = next_char_first_sa_index + sa_start_offset;
  auto new_end = next_char_first_sa_index + sa_end_offset - 1;
  return SA_Interval{new_start, new_end};
}

SearchStates gram::search_base_backwards(const int_Base &pattern_char,
                                         SearchStates const &search_states,
                                         const PRG_Info &prg_info) {
  // Compute the first occurrence of `pattern_char` in the suffix array.
  // Necessary for backward search.
  auto char_alphabet_rank = prg_info.fm_index.char2comp[pattern_char];
  auto char_first_sa_index = prg_info.fm_index.C[char_alphabet_rank];

  SearchStates new_search_states;
  for (auto const &search_state : search_states) {
    auto const new_search_state = search_fm_index_base_backwards(
        pattern_char, char_first_sa_index, search_state, prg_info);
    if (new_search_state)
      new_search_states.push_back(std::move(*new_search_state));
  }
  return new_search_states;
}

std::string gram::serialize_search_state(const SearchState &search_state) {
  std::stringstream ss;
  ss << "****** Search State ******" << std::endl;

  ss << "SA interval: [" << search_state.sa_interval.first << ", "
     << search_state.sa_interval.second << "]";
  ss << std::endl;

  if (not search_state.traversed_path.empty()) {
    ss << "Variant site path [marker, allele id]: " << std::endl;
    for (const auto &variant_site : search_state.traversed_path) {
      auto marker = variant_site.first;

      if (variant_site.second != 0) {
        const auto &allele_id = variant_site.second;
        ss << "[" << marker << ", " << allele_id << "]" << std::endl;
      }
    }
  }
  ss << "****** END Search State ******" << std::endl;
  return ss.str();
}

std::ostream &gram::operator<<(std::ostream &os,
                               const SearchState &search_state) {
  os << serialize_search_state(search_state);
  return os;
}
