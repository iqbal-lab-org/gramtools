"""
Words of caution:
    - read should be above --kmer-size used for indexing
    - watch out for reverse complemented mapping reads which can create multi mapping instances.

 To inspect the integers in the binary PRG String:
 ```hexdump -v -e '1/4 "%d "' filename.bin``` (4-byte integers. beware of endianness)
"""

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
data_dir = base_dir / "tests" / "data" / "prgs"

gramtools._setup_parser()


class IntegrationRunner:
    """
    In charge of calling build, quasimap, and loading the coverage produced in quasimap
    """

    def __init__(self, test_data_dir):
        self.prg_file = str(data_dir / test_data_dir / "prg.bin")
        self.reads_file = str(data_dir / test_data_dir / "reads.fastq")

        self.gram_dir = tempfile.mkdtemp()
        self.run_dir = str(Path(self.gram_dir) / "run")

    def __del__(self):
        shutil.rmtree(self.gram_dir)

    def run_build(self):
        args = gramtools.root_parser.parse_args(
            f"build --gram-dir {self.gram_dir} --prg {self.prg_file} --kmer-size 5".split()
        )
        build.run(args)

    def run_quasimap(self):
        args = gramtools.root_parser.parse_args(
            f"quasimap --gram-dir {self.gram_dir} --run-dir {self.run_dir} --reads {self.reads_file}".split()
        )
        quasimap.run(args)

        ## Make paths to the produced coverage ##
        self.coverage_paths = paths.generate_quasimap_paths(args)

    def run_all(self):
        self.run_build()
        self.run_quasimap()

    def get_pb_cov(self):
        all_per_base_coverage = infer._load_per_base_coverage(
            self.coverage_paths["allele_base_coverage"]
        )
        return all_per_base_coverage

    def get_gped_allele_counts(self):
        allele_groups, all_groups_site_counts = infer._load_grouped_allele_coverage(
            self.coverage_paths["grouped_allele_counts_coverage"]
        )
        return allele_groups, all_groups_site_counts


class twoSites_twoReads_NoNesting(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Setup the attributes only once and share across the test cases

        PRG : "AAA5CC6TA5AC7TTTT8GGG7"
        Reads : "AAATAACGG" and "CACTTTT"
        """
        cls.test_data_dir = "IT1"

        cls.integrator = IntegrationRunner(cls.test_data_dir)

        # Expected per base cov
        cls.expected_per_base = [[[0, 1], [1, 1]], [[1, 1, 1, 1], [1, 1, 0]]]

        # Expect one read mapping to each allele of each of the two sites in the PRG
        cls.expected_grouped_counts = [{"1": 1, "0": 1}, {"1": 1, "0": 1}]

    def test_all(self):
        self.integrator.run_all()

        all_per_base_coverage = self.integrator.get_pb_cov()
        self.assertEqual(self.expected_per_base, all_per_base_coverage)

        allele_groups, all_groups_site_counts = self.integrator.get_gped_allele_counts()
        # Test grouped allele counts
        # First make sure the group IDs are right #
        self.assertEqual(allele_groups["0"], {0})
        self.assertEqual(allele_groups["1"], {1})
        self.assertEqual(self.expected_grouped_counts, all_groups_site_counts)


class twoSites_twoReadswithEquivClasses_NoNesting(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        PRG : "TT5AAAc6AAAg6gg7cAA8gAA8TTCAA"
        Reads : "TTAAA"; "AATTCAA"
        """
        cls.test_data_dir = "IT2"

        cls.integrator = IntegrationRunner(cls.test_data_dir)

        # Expected per base cov
        cls.expected_per_base = [[[1, 1, 1, 0], [1, 1, 1, 0]], [[0, 1, 1], [0, 1, 1]]]

        # Expect each read to map to both alleles of one site each.
        cls.expected_grouped_counts = [{"0": 1}, {"0": 1}]

    def test_all(self):
        self.integrator.run_all()

        all_per_base_coverage = self.integrator.get_pb_cov()
        self.assertEqual(self.expected_per_base, all_per_base_coverage)

        allele_groups, all_groups_site_counts = self.integrator.get_gped_allele_counts()
        self.assertEqual(len(allele_groups), 1)
        self.assertEqual(allele_groups["0"], {0, 1})
        self.assertEqual(self.expected_grouped_counts, all_groups_site_counts)


class OneRead_SNPNestedInsideDeletion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        PRG : "T[cCCC[A,g]CT,]ATTTTt"
        Reads : "CCCAC"; "TATTTT"
        """
        cls.test_data_dir = "IT3"
        # Note that "TATTTT" matches the direct deletion (allele 2) AND
        # inside allele 1, of the outer site.

        cls.integrator = IntegrationRunner(cls.test_data_dir)

        # No expected per base cov: it gets captured in a graph

        # Expect the 'nested' read to match to first allele of both sites,
        # and the 'deletion' read to match the second allele of the first site.
        cls.expected_grouped_counts = [{"0": 1, "1": 1}, {"1": 1}]

    def test_all(self):
        self.integrator.run_all()

        allele_groups, all_groups_site_counts = self.integrator.get_gped_allele_counts()
        self.assertEqual(allele_groups["0"], {0, 1})
        self.assertEqual(allele_groups["1"], {0})
        self.assertEqual(self.expected_grouped_counts, all_groups_site_counts)


if __name__ == "__main__":
    unittest.main(buffer=True)
