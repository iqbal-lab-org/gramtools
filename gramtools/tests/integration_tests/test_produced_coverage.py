"""
Words of caution:
    - read should be above --kmer-size used for indexing
    - watch out for reverse complemented mapping reads which can create multi mapping instances.

 To inspect the integers in the binary PRG String:
 ```hexdump -v -e '1/4 "%d "' filename.bin``` (4-byte integers. beware of endianness)
"""

from gramtools import gramtools_main
from gramtools.commands import paths
from gramtools.commands.build import build
from gramtools.commands.genotype import genotype, utils

import tempfile
import shutil
from pathlib import Path
import unittest

base_dir = Path(__file__ + "/../../../").resolve()
data_dir = base_dir / "tests" / "data" / "prgs"

gramtools_main._setup_parser()


class IntegrationRunner:
    """
    In charge of calling build, quasimap, and loading the coverage produced in quasimap
    """

    def __init__(self, test_data_dir):
        self.prg_file = str(data_dir / test_data_dir / "prg.bin")
        self.reads_file = str(data_dir / test_data_dir / "reads.fastq")

        self.gram_dir = tempfile.mkdtemp()
        self.geno_dir = str(Path(self.gram_dir) / "run")
        self.geno_paths = paths.GenotypePaths(self.geno_dir)

    def __del__(self):
        shutil.rmtree(self.gram_dir)

    def run_build(self):
        args = gramtools_main.root_parser.parse_args(
            f"build --gram_dir {self.gram_dir} --prg {self.prg_file} "
            f" --reference None --kmer_size 5 --force".split()
        )
        build.run(args)

    def run_genotype(self):
        args = gramtools_main.root_parser.parse_args(
            f"genotype --gram_dir {self.gram_dir} --genotype_dir "
            f"{self.geno_dir} --reads {self.reads_file} --sample_id test --force".split()
        )
        genotype.run(args)

    def run_all(self):
        self.run_build()
        self.run_genotype()

    def get_pb_cov(self):
        all_per_base_coverage = utils._load_per_base_coverage(self.geno_paths.pb_cov)
        return all_per_base_coverage

    def get_gped_allele_counts(self):
        allele_groups, all_groups_site_counts = utils._load_grouped_allele_coverage(
            self.geno_paths.gped_cov
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

        # Expect the 'nested' read to match to first allele of both sites,
        # and the 'deletion' read to match the second allele of the first site.
        cls.expected_grouped_counts = [{"0": 1, "1": 1}, {"1": 1}]

    def test_all(self):
        self.integrator.run_all()

        # The PRG has nested sites, so no expected per base coverage recorded
        pb_cov = self.integrator.get_pb_cov()
        self.assertEqual(pb_cov, [])

        allele_groups, all_groups_site_counts = self.integrator.get_gped_allele_counts()
        self.assertEqual(allele_groups["0"], {0, 1})
        self.assertEqual(allele_groups["1"], {0})
        self.assertEqual(self.expected_grouped_counts, all_groups_site_counts)


if __name__ == "__main__":
    unittest.main(buffer=True)
