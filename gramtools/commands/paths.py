## @file
# Sets up file names and directory structures for storing the gramtools inputs and outputs.
import os
import logging
import shutil
from pathlib import Path
from typing import List
from abc import abstractmethod, ABCMeta

log = logging.getLogger("gramtools")


class ProjectPaths(metaclass=ABCMeta):
    @staticmethod
    def path_fact(base_dir: Path):
        """
        Factory building functions which build paths from a base_dir.
       """
        return lambda fname: base_dir / fname

    all_vars = {}

    def check_exists(self, fname: Path, file_description="File"):
        if not fname.exists():
            error_message = f"{file_description} required but not found: {fname}"
            log.error(error_message)
            self.cleanup()
            exit(1)

    @abstractmethod
    def initial_setup(self):
        pass

    @abstractmethod
    def cleanup(self):
        pass

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
    def __init__(self, gram_dir):
        self.gram_dir = Path.resolve(Path(gram_dir))
        self.build_path = ProjectPaths.path_fact(self.gram_dir)

        self.prg = self.build_path("prg")
        self.ref = self.build_path("original_reference.fasta")
        self.built_vcf = self.build_path("build.vcf")
        self.report = self.build_path("build_report.json")
        self.fm_index = self.build_path("fm_index")
        self.cov_graph = self.build_path("cov_graph")
        self.made_gram_dir = False

        self.initial_setup()

    def initial_setup(self):
        if not self.gram_dir.exists():
            log.debug("Creating gram directory:\n%s", self.gram_dir)
            self.gram_dir.mkdir()
            self.made_gram_dir = True

    def cleanup(self):
        if self.gram_dir.exists() and self.made_gram_dir:
            self.gram_dir.rmdir()

    def _check_ref_and_make_ref_link(self, ref_path):
        resolved_ref_path = Path.resolve(Path(ref_path))
        self.check_exists(resolved_ref_path)
        if os.path.lexists(self.ref):
            self.ref.unlink()
        os.symlink(resolved_ref_path, self.ref)

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
        self.force = force  # Whether to erase existing geno_dir

        self.initial_setup()

    def initial_setup(self):
        if self.geno_dir.exists():
            if not self.force:
                log.error(
                    f"{self.geno_dir} already exists.\n"
                    f"Pass --force at command-line to erase it."
                )
                exit(1)
            shutil.rmtree(self.geno_dir)
        self.geno_dir.mkdir()
        self.reads_dir.mkdir()

    def cleanup(self):
        shutil.rmtree(self.geno_dir)

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
            target.symlink_to(read_file)


class SimulatePaths(ProjectPaths):
    def __init__(self, output_dir, sample_id: str, prg_filepath, force=False):
        self.output_dir = Path(output_dir).resolve()
        self.made_output_dir = False
        self.force = force

        self.prg_fpath = Path(prg_filepath).resolve()
        self.json_out = self.output_dir / f"{sample_id}.json"
        self.fasta_out = self.output_dir / f"{sample_id}.fasta"

        self.initial_setup()

    def initial_setup(self):
        self.check_exists(self.prg_fpath)

        if not self.output_dir.exists():
            self.output_dir.mkdir()
            self.made_output_dir = True
        else:
            for path in [self.json_out, self.fasta_out]:
                if path.exists() and not self.force:
                    self.raise_error(
                        f"{path} already exists.\n" f"Run with --force to overwrite."
                    )

    def cleanup(self):
        if self.made_output_dir:
            shutil.rmtree(self.output_dir)


def generate_infer_paths(args):
    all_paths = _generate_run_paths(args.run_dir)
    build_dir = all_paths[
        "build_dir"
    ]  #  Symlink to original --gram-dir, index against which quasimap was ran.
    all_paths.update(_generate_project_paths(build_dir))

    if not os.path.exists(all_paths["infer_dir"]):
        os.mkdir(all_paths["infer_dir"])

    return all_paths


## Make `discover` file paths. Includes making the paths to the reads files used in `quasimap`.
def generate_discover_paths(args):

    all_paths = _generate_run_paths(args.run_dir)
    build_dir = all_paths[
        "build_dir"
    ]  #  Symlink to original --gram-dir, index against which quasimap and infer ran.
    all_paths.update(_generate_project_paths(build_dir))

    #  Make discovery output dir
    if not os.path.exists(all_paths["discover_dir"]):
        os.mkdir(all_paths["discover_dir"])

    #  Check `infer` files are present
    _check_exists(
        all_paths["inferred_fasta"],
        file_description="Fasta formatted inferred personalised reference",
    )

    _check_exists(
        all_paths["inferred_vcf"],
        file_description="Vcf formatted inferred personalised reference",
    )

    #  Build read file paths by scanning the reads_directory, and following symlinks.
    if args.reads is None:
        reads_files = []
        with os.scandir(all_paths["reads_dir"]) as it:
            for entry in it:
                reads_files.append(
                    os.path.realpath(entry.path)
                )  # realpath resolves symlinks, unlike abspath

        all_paths.update({"reads_files": reads_files})

    return all_paths


## Generates paths that will be shared between 'run' commands.
def _generate_run_paths(run_dir):

    # Infer paths
    infer_path = path_fact(paths["infer_dir"])
    paths.update(
        {
            "inferred_fasta": infer_path("inferred.fasta"),
            "inferred_vcf": infer_path("inferred.vcf"),
            "inferred_ref_size": infer_path("inferred_ref_size"),
        }
    )

    # Discover paths
    discover_path = path_fact(paths["discover_dir"])
    paths.update(
        {
            "cortex_vcf": discover_path("cortex.vcf"),
            "rebased_vcf": discover_path("rebased.vcf"),
        }
    )

    return paths
