# To inspect the integers in the binary PRG String:
# ```hexdump -v -e '1/4 "%d "' filename.bin``` (4-byte integers. beware of endianness)

from gramtools import gramtools
from gramtools import paths
from gramtools.commands.build import build
from gramtools.commands import quasimap
from gramtools.commands.infer import infer

import tempfile
import shutil
from pathlib import Path
import unittest

base_dir = Path(__file__).parent.parent.parent
data_dir = base_dir / 'tests' / 'data' / 'prgs'

gramtools._setup_parser()


class twoSitesNoNesting(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Setup the attributes only once and share across the test cases

        The PRG is : "AAA5CC6TA5AC7TTTT8GGG7"
        The reads are : "AAATAACGG" and "CACTTTT"
        """
        cls.gram_dir = tempfile.mkdtemp()

        cls.run_dir = str(Path(cls.gram_dir) / 'run')
        cls.prg_file = str(data_dir / 'two_sites_noNesting.fasta.max_nest5.min_match1.bin')
        cls.reads_file = str(data_dir / 'two_sites_noNesting_twoReads.fastq')

        # Expected per base cov
        cls.expected_per_base = [
            [[0, 1], [1, 1]],
            [[1, 1, 1, 1], [1, 1, 0]]
        ]

        # Expect one read mapping to each allele of each of the two sites in the PRG
        cls.expected_grouped_counts = [
            {'1' : 1, '0' : 1},
            {'1' : 1, '0' : 1},
        ]


    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.gram_dir)

    def test_build(self):
        args = gramtools.root_parser.parse_args(
            f"build --gram-dir {self.gram_dir} --prg {self.prg_file} --kmer-size 5 --all-kmers".split()
        )
        build.run(args)

    def test_quasimap(self):
        args = gramtools.root_parser.parse_args(
            f"quasimap --gram-dir {self.gram_dir} --run-dir {self.run_dir} --reads {self.reads_file}".split()
        )
        quasimap.run(args)

        ## Now test the produced coverage ##
        _paths = paths.generate_quasimap_paths(args)

        # Load coverage stats
        all_per_base_coverage = infer._load_per_base_coverage(_paths['allele_base_coverage'])
        allele_groups, all_groups_site_counts = \
            infer._load_grouped_allele_coverage(_paths['grouped_allele_counts_coverage'])

        self.assertEqual(self.expected_per_base, all_per_base_coverage)

        # Test grouped allele counts
        # First make sure the group IDs are right #
        self.assertEqual(allele_groups['0'], {0})
        self.assertEqual(allele_groups['1'], {1})
        self.assertEqual(self.expected_grouped_counts, all_groups_site_counts)


if __name__ == "__main__":
    unittest.main(buffer=True)
