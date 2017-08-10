import unittest

from . import common
from .. import prg
from .. import kmers


class TestDirectionalRegionRange(unittest.TestCase):
    def _analyse_case(self, start_region_idx, kmer_region_size,
                      prg_structure, expected, reverse):
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)
        start_region = regions[start_region_idx]
        ordered_alleles = kmers._directional_alleles_range(kmer_region_size,
                                                           start_region,
                                                           regions,
                                                           reverse=reverse)
        result = list(ordered_alleles)
        self.assertEqual(result, expected)

    def test_rightEdgeStartingRegion_noForwardRegions(self):
        prg_structure = [
            ['C'],
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
        ]
        expected = []
        start_region_idx = -1
        kmer_region_size = 100
        reverse = False
        self._analyse_case(start_region_idx, kmer_region_size,
                           prg_structure, expected, reverse)

    def test_leftEdgeStartingRegion_noReverseRegions(self):
        prg_structure = [
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
            ['C'],
        ]
        expected = []
        start_region_idx = 0
        kmer_region_size = 100
        reverse = True
        self._analyse_case(start_region_idx, kmer_region_size,
                           prg_structure, expected, reverse)

    def test_populatedRevRegion_revRegionCorrectlyOrdered(self):
        prg_structure = [
            ['C', 'G'],
            ['G'],
            ['A', 'T'],
            ['CCCCCCCC'],
            ['GC', 'TC'],
        ]
        expected = [
            ['C', 'G'],
            ['G'],
        ]
        start_region_idx = 2
        kmer_region_size = 2
        reverse = True
        self._analyse_case(start_region_idx, kmer_region_size,
                           prg_structure, expected, reverse)

    def test_populatedFwdRegion_fwdRegionCorrectlyOrdered(self):
        prg_structure = [
            ['G'],
            ['GC', 'TC'],
            ['C'],
            ['A', 'T'],
            ['G'],
            ['C', 'G'],
        ]
        expected = [
            ['G'],
            ['C', 'G'],
        ]
        start_region_idx = 3
        kmer_region_size = 2
        reverse = False
        self._analyse_case(start_region_idx, kmer_region_size,
                           prg_structure, expected, reverse)


class TestRegionsWithinDistance(unittest.TestCase):
    def _analyse_case(self, start_region_idx, kmer_region_size,
                      prg_structure, expected):
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)
        start_region = regions[start_region_idx]
        ordered_alleles = kmers._alleles_within_range(kmer_region_size,
                                                      start_region,
                                                      regions)
        result = list(ordered_alleles)
        self.assertEqual(result, expected)

    def test_distanceExcludesRegion_lastRegionExcluded(self):
        prg_structure = [
            ['A', 'C'],
            ['CC'],
            ['G', 'T'],
        ]
        expected = [
            ['A', 'C'],
            ['C'],
        ]
        start_region_idx = 0
        kmer_region_size = 1
        self._analyse_case(start_region_idx,
                           kmer_region_size,
                           prg_structure,
                           expected)

    def test_thresholdMaxDistance_startingRegionDistanceIgnored(self):
        prg_structure = [
            ['AAAAAAA', 'CCCCC'],
            ['CC'],
            ['G', 'T'],
        ]
        expected = [
            ['AAAAAAA', 'CCCCC'],
            ['CC'],
        ]
        start_region_idx = 0
        kmer_region_size = 2
        self._analyse_case(start_region_idx,
                           kmer_region_size,
                           prg_structure,
                           expected)

    def test_distanceEdgeIncludeLargeRegion_largeRegionInRange(self):
        prg_structure = [
            ['G'],
            ['A', 'T'],
            ['CC'],
            ['GGGGGGGG', 'T'],
        ]
        expected = [
            ['G'],
            ['A', 'T'],
            ['CC'],
            ['G', 'T'],
        ]
        start_region_idx = 1
        kmer_region_size = 3
        self._analyse_case(start_region_idx,
                           kmer_region_size,
                           prg_structure,
                           expected)

    def test_shortDistance_excludesFwdRevRegions(self):
        prg_structure = [
            ['C', 'G'],
            ['G'],
            ['A', 'T'],
            ['CC'],
            ['GC', 'TC'],
        ]
        expected = [
            ['G'],
            ['A', 'T'],
            ['C'],
        ]
        start_region_idx = 2
        kmer_region_size = 1
        self._analyse_case(start_region_idx,
                           kmer_region_size,
                           prg_structure,
                           expected)

    def test_largeFwdRegion_postFwdRegionsExcluded(self):
        prg_structure = [
            ['C', 'G'],
            ['G'],
            ['A', 'T'],
            ['CCCCCCCC'],
            ['GC', 'TC'],
        ]
        expected = [
            ['C', 'G'],
            ['G'],
            ['A', 'T'],
            ['CC'],
        ]
        start_region_idx = 2
        kmer_region_size = 2
        self._analyse_case(start_region_idx,
                           kmer_region_size,
                           prg_structure,
                           expected)

    def test_largeRevRegion_preRevRegionsExcluded(self):
        prg_structure = [
            ['T', 'A'],
            ['CCCCCCCC'],
            ['C', 'G'],
            ['G'],
            ['A', 'T'],
            ['GGGGGGG'],
        ]
        expected = [
            ['C', 'G'],
            ['G'],
            ['A', 'T'],
            ['GG'],
        ]
        start_region_idx = 4
        kmer_region_size = 2
        self._analyse_case(start_region_idx,
                           kmer_region_size,
                           prg_structure,
                           expected)


class TestGenomePaths(unittest.TestCase):
    def _analyse_case(self, start_region_idx, kmer_region_size,
                      prg_structure, expected):
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)
        start_region = regions[start_region_idx]
        ordered_alleles = kmers._alleles_within_range(kmer_region_size,
                                                      start_region,
                                                      regions)
        results = [path for path in
                   kmers._genome_paths(ordered_alleles)]
        self.assertEqual(results, expected)

    def test_twoVariantRegionsInRange_pathsIncludeAllAlleles(self):
        prg_structure = [
            ['AC'],
            ['T', 'G'],
            ['C'],
            ['A', 'T'],
        ]
        expected = [
            'ACTCA',
            'ACTCT',
            'ACGCA',
            'ACGCT',
        ]
        start_region_idx = 1
        kmer_region_size = 3
        self._analyse_case(start_region_idx, kmer_region_size,
                           prg_structure, expected)

    def test_largeEdgeAllele_alleleTruncated(self):
        prg_structure = [
            ['AC'],
            ['T', 'G'],
            ['C'],
            ['AGGGG', 'T'],
        ]
        expected = [
            'ACTCA',
            'ACTCT',
            'ACGCA',
            'ACGCT',
        ]
        start_region_idx = 1
        kmer_region_size = 2
        self._analyse_case(start_region_idx, kmer_region_size,
                           prg_structure, expected)

    def test_largeKmerRegionSize_correctPathsGenerated(self):
        prg_structure = [
            ['AC'],
            ['T', 'G'],
            ['CT'],
        ]
        expected = [
            'ACTCT',
            'ACGCT',
        ]
        start_region_idx = 1
        kmer_region_size = 1000
        self._analyse_case(start_region_idx, kmer_region_size,
                           prg_structure, expected)

    def test_manyAllelesInVarRegion_pathsIncludeAllAlleles(self):
        prg_structure = [
            ['AC'],
            ['T', 'G', 'C', 'AAA'],
            ['GT'],
        ]
        expected = [
            'ACTGT',
            'ACGGT',
            'ACCGT',
            'ACAAAGT',
        ]
        start_region_idx = 1
        kmer_region_size = 2
        self._analyse_case(start_region_idx, kmer_region_size,
                           prg_structure, expected)


class TestKmersFromGenomePaths(unittest.TestCase):
    def _analyse_case(self, genome_paths, kmer_size, expected):
        results = [kmer for kmer in
                   kmers._kmers_from_genome_paths(genome_paths, kmer_size)]
        self.assertEqual(results, expected)

    def test_singleGenomePath_allKmersGenerated(self):
        genome_paths = [
            'ACTGT',
        ]
        expected = [
            'ACT',
            'CTG',
            'TGT',
        ]
        kmer_size = 3
        self._analyse_case(genome_paths, kmer_size, expected)

    def test_multipleGenomePaths_allKmersGenerated(self):
        genome_paths = [
            'ACTGT',
            'TATC',
        ]
        expected = [
            'ACT',
            'CTG',
            'TGT',
            'TAT',
            'ATC',
        ]
        kmer_size = 3
        self._analyse_case(genome_paths, kmer_size, expected)

    def test_repeatedBases_duplicateKmersExcluded(self):
        genome_paths = [
            'TTTTTTTTTTT',
        ]
        expected = [
            'TTT',
        ]
        kmer_size = 3
        self._analyse_case(genome_paths, kmer_size, expected)


class TestGenerate(unittest.TestCase):
    def _analyse_case(self, expected, prg_structure, kmer_size,
                      kmer_region_size, nonvariant_kmers):
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)
        generated_kmers = kmers._generate(kmer_region_size, kmer_size,
                                          regions, nonvariant_kmers)
        result = list(generated_kmers)
        self.assertEqual(result, expected)

    def test_complexPrgStructure_allKmersGenerated(self):
        prg_structure = [
            ['TC', 'A'],
            ['CTTG'],
        ]
        kmer_size = 3
        kmer_region_size = 2
        nonvariant_kmers = False

        expected = [
            'TCC',
            'CCT',
            'ACT',
        ]
        self._analyse_case(expected, prg_structure, kmer_size,
                           kmer_region_size, nonvariant_kmers)

    def test_allelesLargerThankmerRegionBoundary_allelesTruncated(self):
        prg_structure = [
            ['AGGCTT'],
            ['A', 'C'],
            ['AGA'],
        ]
        kmer_size = 3
        kmer_region_size = 2
        nonvariant_kmers = False

        expected = [
            'TTA',
            'TAA',
            'AAG',
            'TTC',
            'TCA',
            'CAG',
        ]
        self._analyse_case(expected, prg_structure, kmer_size,
                           kmer_region_size, nonvariant_kmers)

    def test_nonvariantKmersGenerated_allKmersGenerated(self):
        prg_structure = [
            ['AGGCTT'],
            ['A', 'C'],
            ['AGA'],
        ]
        kmer_size = 3
        kmer_region_size = 2
        nonvariant_kmers = True

        expected = [
            'AGG',
            'GGC',
            'GCT',
            'CTT',
            'TTA',
            'TAA',
            'TTC',
            'TCA',
            'AAG',
            'CAG',
            'AGA',
        ]
        self._analyse_case(expected, prg_structure, kmer_size,
                           kmer_region_size, nonvariant_kmers)
