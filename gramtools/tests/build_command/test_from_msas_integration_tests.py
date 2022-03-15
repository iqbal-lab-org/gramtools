from unittest import TestCase
from pathlib import Path
from tempfile import mkdtemp
from shutil import rmtree
from itertools import chain

from pybedtools import BedTool, Interval

from gramtools.tests import data_dir
from gramtools.commands.build.from_msas import standalone_build_from_msas
from gramtools.commands.common import load_fasta, ints_to_bytes, nuc_translation

data_dir = data_dir / "from_msas"


class SetupFiles:
    def __init__(self, base_dir: Path):
        # File names in to_build.bed are relative to where to_build.bed is, for portability of test code.
        # Here we convert them into absolute file names so they can actually be found
        self.built_prg_dirname = mkdtemp()
        intervals = BedTool(str(base_dir / "to_build.bed"))
        new_intervals = list()
        for interval in intervals:
            new_name = str((base_dir / interval.name).resolve())
            new_intervals.append(
                Interval(interval.chrom, interval.start, interval.end, new_name)
            )
        self.prgs_bed = f"{self.built_prg_dirname}/to_build_absolute_fnames.bed"
        BedTool(new_intervals).saveas(self.prgs_bed)
        self.reference = base_dir / "ref.fa"
        self.coords_file = base_dir / "chrom_sizes.tsv"

    def __del__(self):
        rmtree(self.built_prg_dirname)


class TestBuildFromPrgsBed(TestCase):
    def test_build_prg_from_prg_bed(self):
        input_files = SetupFiles(data_dir)
        chrom_seqs = load_fasta(input_files.reference)
        all_ref_nucleotides = list(chain.from_iterable(map(list, chrom_seqs.values())))
        expected_prg_ints = [nuc_translation[nuc] for nuc in all_ref_nucleotides]

        # Add in the variant I created in the test integration data (CC vs CA in second chromosome)
        expected_prg_ints = (
            expected_prg_ints[:6]
            + [
                5,
                nuc_translation["C"],
                nuc_translation["C"],
                6,
                nuc_translation["C"],
                nuc_translation["A"],
                6,
            ]
            + expected_prg_ints[8:]
        )
        expected_prg = ints_to_bytes(expected_prg_ints)
        result_built_bed, result_prg_ints = standalone_build_from_msas(
            input_files.prgs_bed,
            input_files.reference,
            input_files.coords_file,
            input_files.built_prg_dirname,
        )
        self.assertEqual(expected_prg_ints, result_prg_ints)

        # Check we built the right bed output. Important to validate is the start/end coords
        expected_start_ends = [(0, 2), (2, 4), (4, 5), (0, 1), (1, 3), (3, 5)]
        actual_start_ends = list(
            map(lambda interval: (interval.start, interval.end), result_built_bed)
        )
        self.assertEqual(expected_start_ends, actual_start_ends)
