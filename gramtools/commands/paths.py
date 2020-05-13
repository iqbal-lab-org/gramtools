## @file
# Sets up file names and directory structures for storing the gramtools inputs and outputs.
import os
import logging
import shutil
from pathlib import Path
from typing import List

log = logging.getLogger("gramtools")


class ProjectPaths:
    @staticmethod
    def path_fact(base_dir: Path):
        """
        Factory building functions which build paths from a base_dir.
       """
        return lambda fname: base_dir / fname

    def __init__(self, output_dir, force: bool):
        self.output_dir = output_dir
        self.made_output_dir = False
        self.force = force  # Whether to erase existing geno_dir

    all_vars = {}

    def check_exists(self, fname: Path, file_description="File"):
        if not fname.exists():
            error_message = f"{file_description} required but not found: {fname}"
            log.error(error_message)
            self.cleanup()
            exit(1)

    def initial_setup(self):
        if not self.output_dir.exists():
            log.debug("Creating output directory:\n%s", self.output_dir)
            self.output_dir.mkdir(parents=True)
            self.made_output_dir = True
            return

        if not self.force:
            self.raise_error(
                f"{self.output_dir} already exists.\n" f"Run with --force to overwrite."
            )
        shutil.rmtree(self.output_dir)
        self.output_dir.mkdir()

    def cleanup(self):
        if self.made_output_dir:
            shutil.rmtree(self.output_dir)

    def isOutputtablePath(self, object):
        if isinstance(object, list):
            for el in object:
                if not isinstance(el, Path):
                    return False
        if not isinstance(object, Path):
            return False

        return True

    def raise_error(self, err_message):
        self.cleanup()
        log.error(err_message)
        exit(1)

    def items(self):
        all_vars = vars(self)
        ProjectPaths.all_vars = {
            var: entry
            for var, entry in all_vars.items()
            if self.isOutputtablePath(entry)
        }
        return ProjectPaths.all_vars.items()

    def dict(self):
        """
        Like items() method but:
        i)returns a dictionary, not its items
        ii)the item Path entries are all stringified
        """
        all_vars = vars(self)
        returned_vars = {}
        for var, entry in all_vars.items():
            if not self.isOutputtablePath(entry):
                continue
            if isinstance(entry, list):
                entry = [str(element) for element in entry]
            else:
                entry = str(entry)
            returned_vars[var] = entry
        return returned_vars


class BuildPaths(ProjectPaths):
    def __init__(self, gram_dir, force=False):
        self.gram_dir = Path.resolve(Path(gram_dir))
        super().__init__(self.gram_dir, force)

        self.build_path = ProjectPaths.path_fact(self.gram_dir)
        self.prg = self.build_path("prg")
        self.coords_file = self.build_path("prg_coords.tsv")
        self.built_vcf = self.build_path("build.vcf")
        self.report = self.build_path("build_report.json")
        self.fm_index = self.build_path("fm_index")
        self.cov_graph = self.build_path("cov_graph")

    def setup(self):
        super().initial_setup()

    def _check_ref_and_make_ref_link(self, ref_path):
        resolved_ref_path = Path.resolve(Path(ref_path))
        self.check_exists(resolved_ref_path)
        self.ref = resolved_ref_path

    def _check_and_flatten_vcf_filenames(self, vcf: List[List[str]]):
        vcf_files = [Path(vcf_file) for arglist in vcf for vcf_file in arglist]

        for vcf_file in vcf_files:
            self.check_exists(vcf_file)
        self.input_vcfs = vcf_files

    def ready_ref_and_vcf(self, reference: str, vcfs: List[List[str]]):
        self._check_ref_and_make_ref_link(reference)
        self._check_and_flatten_vcf_filenames(vcfs)


class GenotypePaths(ProjectPaths):
    def __init__(self, genotype_dir, force=False):
        self.geno_dir = Path.resolve(Path(genotype_dir))
        super().__init__(self.geno_dir, force)

        self.geno_path = ProjectPaths.path_fact(self.geno_dir)
        self.cov_path = ProjectPaths.path_fact(self.geno_dir / "coverage")
        self.results_path = ProjectPaths.path_fact(self.geno_dir / "genotype")

        self.gram_dir = self.geno_path("gram_dir")
        self.reads_dir = self.geno_path("reads_dir")
        self.report = self.geno_path("genotype_report.json")

        # Quasimapping-related
        self.read_stats = self.geno_path("read_stats.json")
        self.gped_cov = self.cov_path("grouped_allele_counts_coverage.json")
        self.pb_cov = self.cov_path("allele_base_coverage.json")

        # Infer-related
        self.geno_vcf = self.results_path("genotyped.vcf.gz")
        self.pers_ref = self.results_path("personalised_reference.fasta")

    def setup(self, args):
        super().initial_setup()
        self.reads_dir.mkdir()
        self._link_to_build(args.gram_dir)
        self._link_to_reads(args.reads)

    def _link_to_build(self, existing_gram_dir):
        """
        Make a reference to the gram_dir made in build. This avoids downstream commands
        to require a --gram_dir argument.
        """
        target = Path(existing_gram_dir).resolve()
        self.check_exists(target)
        if os.path.lexists(self.gram_dir):
            os.unlink(self.gram_dir)
        self.gram_dir.symlink_to(target, target_is_directory=True)

    def _link_to_reads(self, reads: List[List[str]]):
        self.reads_files = [
            Path(read_file).resolve() for arglist in reads for read_file in arglist
        ]

        #  Make symlinks to the read files
        for read_file in self.reads_files:
            target = (
                self.reads_dir / read_file.name
            )  # name attrib is the file's os.basename
            self.check_exists(read_file)
            target.symlink_to(read_file)


class DiscoverPaths(ProjectPaths):
    def __init__(self, discovery_dir, genotype_dir, force=False):
        self.disco_dir = Path(discovery_dir).resolve()
        super().__init__(self.disco_dir, force)

        # Build paths from genotype dir
        geno_paths = GenotypePaths(genotype_dir)
        self.pers_ref = geno_paths.pers_ref
        self.geno_vcf = geno_paths.geno_vcf
        self.geno_report = geno_paths.report

        self.reads_files = []
        self.check_exists(geno_paths.reads_dir)
        for read_file in geno_paths.reads_dir.iterdir():
            self.reads_files.append(
                read_file.resolve()
            )  # to absolute path + resolves symlinks

        self.discov_vcf_cortex = self.disco_dir / "cortex.vcf"
        self.final_vcf = self.disco_dir / "final.vcf"
        self.rebasing_map = self.disco_dir / "rebasing_map.json"

    def setup(self):
        super().initial_setup()
        self.check_exists(self.pers_ref)


class SimulatePaths(ProjectPaths):
    def __init__(self, output_dir, sample_id: str, prg_filepath, force=False):
        self.sim_dir = Path(output_dir).resolve()
        super().__init__(self.sim_dir, force)

        self.prg_fpath = Path(prg_filepath).resolve()
        self.json_out = self.sim_dir / f"{sample_id}.json"
        self.fasta_out = self.sim_dir / f"{sample_id}.fasta"

    def setup(self):
        if not self.sim_dir.exists():
            self.sim_dir.mkdir()
            self.made_output_dir = True

        self.check_exists(self.prg_fpath)
        for path in [self.json_out, self.fasta_out]:
            if path.exists() and not self.force:
                self.raise_error(
                    f"{self.sim_dir} already exists.\n"
                    f"Run with --force to overwrite."
                )
