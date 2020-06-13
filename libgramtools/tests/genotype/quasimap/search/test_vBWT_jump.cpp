/**
 * @file
 * Unit tests for vBWT backward searching.
 * Terminology:
 *  - A variant locus is where you find variant **markers**;
 *  = pairs of site & allele markers.
 *  - A site 'entry' (resp. 'exit') is the 3' (resp. 5') part
 *  of a site in the linear PRG; because we are mapping backwards.
 *
 * Test suites:
 *  - MarkerSearch: checking finding and positioning variant markers in the PRG
 * string
 *  - MarkerSAIntervals: Recovering SA Interval of variant markers.
 *  - VariantLocus_Path: checking search recovers right variant site/allele
 * combinations.
 *
 *  - SearchStateJump: vBWT jumping producing correct new `SearchState`s
 *  - SearchStateJump_Nested: same as above, on nested PRG strings
 */
#include <cctype>

#include "gtest/gtest.h"

#include "build/kmer_index/build.hpp"
#include "genotype/quasimap/search/vBWT_jump.hpp"
#include "prg/prg_info.hpp"
#include "submod_resources.hpp"

using namespace gram;

/*
PRG: GCGCT5C6G6A6AGTCCT
i	BWT	SA	text_suffix
0	G	18
1	6	12	A G T C C T
2	6	10	A 6 A G T C C T
3	G	15	C C T
4	T	1	C G C T 5 C 6 G 6 A 6 A G T C C T
5	C	16	C T
6	T	3	C T 5 C 6 G 6 A 6 A G T C C T
7	5	6	C 6 G 6 A 6 A G T C C T
8	0	0	G C G C T 5 C 6 G 6 A 6 A G T C C T
9	C	2	G C T 5 C 6 G 6 A 6 A G T C C T
10	A	13	G T C C T
11	6	8	G 6 A 6 A G T C C T
12	C	17	T
13	T	14	T C C T
14	C	4	T 5 C 6 G 6 A 6 A G T C C T
15	G	5	5 C 6 G 6 A 6 A G T C C T
16	A	11	6 A G T C C T
17	T	9	6 A 6 A G T C C T
18	C	7	6 G 6 A 6 A G T C C T
*/

TEST(MarkerSearch, GivenCharA_FindLeftMarkers_AndSeedSearchStates) {
  auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
  auto prg_info = generate_prg_info(prg_raw);
  // first char: a
  SearchState initial_search_state = {SA_Interval{1, 2}};

  auto result = left_markers_search(initial_search_state, prg_info);
  MarkersSearchResults expected = {
      {6, ALLELE_UNKNOWN},
      {5, FIRST_ALLELE + 2},
  };
  EXPECT_EQ(result, expected);

  // Expect three: one for exiting the site; two for entering.
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);
  EXPECT_EQ(markers_search_states.size(), 2);
}

// The convention is as follows: if the position marks a site exit, the marker
// will be a site marker, and if it marks a site entry, the marker will be an
// allele marker.
TEST(MarkerSearch, TestSiteMarker_Entry_or_Exit) {
  auto prg_raw = encode_prg("gcgct5C6g6a6Agtcct");
  auto prg_info = generate_prg_info(prg_raw);

  // TEST 1: char a at site entry point
  SearchState search_state = {SA_Interval{1, 1}};

  auto result = left_markers_search(search_state, prg_info);

  auto variant_marker = result[0].first;
  EXPECT_TRUE(is_allele_marker(variant_marker));

  // TEST 2: char c at site exit point
  search_state = {SA_Interval{7, 7}};
  result = left_markers_search(search_state, prg_info);
  variant_marker = result[0].first;
  EXPECT_TRUE(is_site_marker(variant_marker));
}

TEST(MarkerSearch, GivenCharG_ReturnOneCorrectSearchResults) {
  auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
  auto prg_info = generate_prg_info(prg_raw);
  // first char: g
  SearchState initial_search_state = {SA_Interval{8, 11}};

  auto result = left_markers_search(initial_search_state, prg_info);
  MarkersSearchResults expected = {
      {5, FIRST_ALLELE + 1},
  };
  EXPECT_EQ(result, expected);
}

TEST(SearchStateJump, SingleCharAllele_CorrectSkipToSiteStartBoundaryMarker) {
  auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
  auto prg_info = generate_prg_info(prg_raw);
  // first char: g
  SearchState initial_search_state = {SA_Interval{8, 11}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);
  const auto &first_markers_search_state = markers_search_states.front();

  const auto &result = first_markers_search_state.sa_interval;
  SA_Interval expected = {15, 15};
  EXPECT_EQ(result, expected);
}

TEST(MarkerSearch, GivenCharG_NoMarkersToLeft) {
  auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
  auto prg_info = generate_prg_info(prg_raw);
  // first char: g
  SearchState initial_search_state = {SA_Interval{8, 11}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);
  const auto &result = markers_search_states.size();
  auto expected = 1;
  EXPECT_EQ(result, expected);
}

TEST(MarkerSearch, GivenCharC_JumptoSiteStart) {
  auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
  auto prg_info = generate_prg_info(prg_raw);
  // first char: c
  SearchState initial_search_state = {SA_Interval{3, 7}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);
  const auto &first_markers_search_state = markers_search_states.front();

  EXPECT_EQ(markers_search_states.size(), 1);
  auto &result = first_markers_search_state.sa_interval;
  SA_Interval expected = {15, 15};
  EXPECT_EQ(result, expected);
}

TEST(MarkerSAIntervals, AlleleMarkerAnd3Alleles_correctSAInterval) {
  auto prg_raw = encode_prg("gcgct5c6g6a6agtcct");
  auto prg_info = generate_prg_info(prg_raw);
  Marker allele_marker = 6;

  auto result = get_allele_marker_sa_interval(allele_marker, prg_info);
  SA_Interval expected = {16, 18};
  EXPECT_EQ(result, expected);
}

TEST(MarkerSAIntervals, AlleleMarkerAnd2Alleles_correctSAInterval) {
  auto prg_raw = encode_prg("aca5g6t6catt");
  auto prg_info = generate_prg_info(prg_raw);

  auto result = get_allele_marker_sa_interval(6, prg_info);
  SA_Interval expected = {11, 12};
  EXPECT_EQ(result, expected);
}

/*
PRG: 7G8C8G9T10A10
i	BWT	SA	text_suffix
0	10	11	1
1	10	9	0 A 1
2	8	3	C 8 G 9 T 1 0 A 1
3	7	1	G 8 C 8 G 9 T 1 0 A 1
4	8	5	G 9 T 1 0 A 1
5	9	7	T 1 0 A 1
6	0	0	7 G 8 C 8 G 9 T 1 0 A 1
7	G	2	8 C 8 G 9 T 1 0 A 1
8	C	4	8 G 9 T 1 0 A 1
9	G	6	9 T 1 0 A 1
10	A	10	A 1
11	T	8	1 0 A 1
*/
TEST(MarkerSAIntervals,
     GivenPrgWithNonContinuousAlphabet_CorrectAlleleMarkerEndBoundary) {
  auto prg_raw = encode_prg("7g8c8g9t10a10");
  auto prg_info = generate_prg_info(prg_raw);

  auto result = get_allele_marker_sa_interval(8, prg_info);
  SA_Interval expected = {7, 8};
  EXPECT_EQ(result, expected);
}

/*
PRG: GCGCT5C6G6T6AGTCCT
i	BWT	SA	text_suffix
0	T	18
1	6	12	A G T C C T
2	T	15	C C T
3	G	1	C G C T 5 C 6 G 6 T 6 A G T C C T
4	C	16	C T
5	G	3	C T 5 C 6 G 6 T 6 A G T C C T
6	5	6	C 6 G 6 T 6 A G T C C T
7	0	0	G C G C T 5 C 6 G 6 T 6 A G T C C T
8	C	2	G C T 5 C 6 G 6 T 6 A G T C C T
9	A	13	G T C C T
10	6	8	G 6 T 6 A G T C C T
11	C	17	T
12	G	14	T C C T
13	C	4	T 5 C 6 G 6 T 6 A G T C C T
14	6	10	T 6 A G T C C T
15	T	5	5 C 6 G 6 T 6 A G T C C T
16	T	11	6 A G T C C T
17	C	7	6 G 6 T 6 A G T C C T
18	G	9	6 T 6 A G T C C T

 */

TEST(SearchStateJump, AtSiteEntry_CorrectSearchStateJump) {
  auto prg_raw = encode_prg("gcgct5c6g6t6Agtcct");
  auto prg_info = generate_prg_info(prg_raw);

  // first char: a
  SearchState initial_search_state = {SA_Interval{1, 1}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);

  EXPECT_EQ(markers_search_states.size(), 1);

  SearchStates expected = {SearchState{
      SA_Interval{16, 18},
      VariantSitePath{},
      VariantSitePath{VariantLocus{5, ALLELE_UNKNOWN}},
  }};

  EXPECT_EQ(markers_search_states, expected);
}

TEST(SearchStateJump, Allele2SiteExit_CorrectSearchStateJump) {
  auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
  auto prg_info = generate_prg_info(prg_raw);

  // first char: g
  SearchState initial_search_state = {SA_Interval{7, 10}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);
  SearchStates expected = {SearchState{
      SA_Interval{15, 15},
      VariantSitePath{VariantLocus{5, FIRST_ALLELE + 1}},
      VariantSitePath{},
  }};
  EXPECT_EQ(markers_search_states, expected);
}

TEST(SearchStateJump, Allele1SiteExit_CorrectSearchStateJump) {
  auto prg_raw = encode_prg("gcgct5c6g6t6agtcct");
  auto prg_info = generate_prg_info(prg_raw);

  // first char: c
  SearchState initial_search_state = {SA_Interval{2, 6}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);
  SearchStates expected = {SearchState{
      SA_Interval{15, 15},
      VariantSitePath{VariantLocus{5, FIRST_ALLELE}},
      VariantSitePath{},
  }};
  EXPECT_EQ(markers_search_states, expected);
}

/**********************/
/* Nested PRG Strings */
/**********************/
/*
PRG: [AC,[C,G]]T
i	BWT	SA	text_suffix
0	T	11	0
1	5	1	A C 6 7 C 8 G 8 6 T 0
2	A	2	C 6 7 C 8 G 8 6 T 0
3	7	5	C 8 G 8 6 T 0
4	8	7	G 8 6 T 0
5	6	10	T 0
6	0	0	5 A C 6 7 C 8 G 8 6 T 0
7	8	9	6 T 0
8	C	3	6 7 C 8 G 8 6 T 0
9	6	4	7 C 8 G 8 6 T 0
10	C	6	8 G 8 6 T 0
11	G	8	8 6 T 0
*/

TEST(SearchStateJump_Nested, DoubleExit_CorrectSearchStateJump) {
  auto prg = prg_string_to_ints("[AC,[C,G]]T");
  auto prg_info = generate_prg_info(prg);

  // first char: c at index 5 in PRG
  SearchState initial_search_state = {SA_Interval{3, 3}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);

  SearchStates expected = {SearchState{
      SA_Interval{6, 6},
      VariantSitePath{VariantLocus{7, FIRST_ALLELE},
                      VariantLocus{5, FIRST_ALLELE + 1}},
      VariantSitePath{},
  }};
  EXPECT_EQ(markers_search_states, expected);
}

TEST(SearchStateJump_Nested, DoubleEntry_CorrectSearchStateJump) {
  auto prg = prg_string_to_ints("[AC,[C,G]]T");
  auto prg_info = generate_prg_info(prg);

  // first char: t
  SearchState initial_search_state = {SA_Interval{5, 5}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);

  SearchStates expected = {SearchState{
                               SA_Interval{7, 8},
                               VariantSitePath{},
                               VariantSitePath{VariantLocus{5, ALLELE_UNKNOWN}},
                           },
                           SearchState{
                               SA_Interval{10, 11},
                               VariantSitePath{},
                               VariantSitePath{
                                   VariantLocus{5, ALLELE_UNKNOWN},
                                   VariantLocus{7, ALLELE_UNKNOWN},
                               },
                           }};
  EXPECT_EQ(markers_search_states, expected);
}

/*
PRG: [C,G][C,G]
i	BWT	SA	text_suffix
0	8	10	0
1	5	1	C 6 G 6 7 C 8 G 8 0
2	7	6	C 8 G 8 0
3	6	3	G 6 7 C 8 G 8 0
4	8	8	G 8 0
5	0	0	5 C 6 G 6 7 C 8 G 8 0
6	C	2	6 G 6 7 C 8 G 8 0
7	G	4	6 7 C 8 G 8 0
8	6	5	7 C 8 G 8 0
9	G	9	8 0
10	C	7	8 G 8 0
 */
TEST(SearchStateJump_Nested, ExitToEntry_CorrectSearchStateJump) {
  auto prg = prg_string_to_ints("[C,G][C,G]");
  auto prg_info = generate_prg_info(prg);

  // first char: c at index 6 in PRG
  SearchState initial_search_state = {SA_Interval{2, 2}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);

  SearchStates expected = {SearchState{
      SA_Interval{6, 7},
      VariantSitePath{VariantLocus{7, FIRST_ALLELE}},
      VariantSitePath{VariantLocus{5, ALLELE_UNKNOWN}},
  }};
  EXPECT_EQ(markers_search_states, expected);
}

/*
PRG: A[C,,G]T
i	BWT	SA	text_suffix
0	T	8	0
1	0	0	A 5 C 6 6 G 6 T 0
2	5	2	C 6 6 G 6 T 0
3	6	5	G 6 T 0
4	6	7	T 0
5	A	1	5 C 6 6 G 6 T 0
6	6	4	6 G 6 T 0
7	G	6	6 T 0
8	C	3	6 6 G 6 T 0
*/

TEST(SearchStateJump_Nested, DirectDeletion_CorrectSearchStateJump) {
  auto prg = prg_string_to_ints("A[C,,G]T");
  auto prg_info = generate_prg_info(prg);

  // first char: T
  // We expect to skip past the direct deletion
  SearchState initial_search_state = {SA_Interval{4, 4}};
  const auto &markers_search_states =
      search_state_vBWT_jumps(initial_search_state, prg_info);

  SearchStates expected = {
      SearchState{
          SA_Interval{6, 8},
          VariantSitePath{},
          VariantSitePath{VariantLocus{5, ALLELE_UNKNOWN}},
      },

      SearchState{
          SA_Interval{5, 5},
          VariantSitePath{VariantLocus{5, FIRST_ALLELE + 1}},
          VariantSitePath{},
      }};

  EXPECT_EQ(markers_search_states, expected);
}
