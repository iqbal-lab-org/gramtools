from tempfile import mkdtemp
from pathlib import Path
import shutil
import json
from unittest import TestCase

from gramtools.commands.simulate import simulate
import gramtools.tests.integration_tests as it_tests

data_dir = it_tests.data_dir


def filter_desc(simu_json):
    for sample_dict in simu_json["Samples"]:
        sample_dict.pop("Desc")


class SimulateRunner:
    def __init__(self, test_data_dir: Path):
        self.prg_file = str(data_dir / test_data_dir / "prg.bin")
        self.simu_dir = Path(mkdtemp())
        self.made_fasta = self.simu_dir / "made.fasta"
        self.made_json = self.simu_dir / "made.json"
        self.induced_json = self.simu_dir / "induced.json"

    def __del__(self):
        shutil.rmtree(self.simu_dir)

    def run_make_paths(self):
        args = it_tests.gramtools_main.root_parser.parse_args(
            f"simulate --prg {self.prg_file} -n 5 --sample_id made -o {self.simu_dir}".split()
        )
        simulate.run(args)

    def run_induce_genotypes(self):
        args = it_tests.gramtools_main.root_parser.parse_args(
            f"simulate --prg {self.prg_file} -o {self.simu_dir} "
            f"--sample_id induced -i {self.made_fasta}".split()
        )
        simulate.run(args)


class TestSimulatedPaths(TestCase):
    def test_make_paths_and_induce_from_paths_get_same_jsons(self):
        runner = SimulateRunner("IT1")
        runner.run_make_paths()
        with runner.made_json.open() as input:
            made_json = json.load(input)
            filter_desc(made_json)

        runner.run_induce_genotypes()
        with runner.induced_json.open() as input:
            induced_json = json.load(input)
            filter_desc(induced_json)

        self.assertEqual(made_json, induced_json)
