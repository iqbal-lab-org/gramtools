## @file
# Sets up file names and directory structures for storing the gramtools inputs and outputs.
import os
import logging
import shutil
from pathlib import Path
from typing import List

log = logging.getLogger("gramtools")


def path_fact(base_dir):
    """
    Factory building functions which build paths from a base_dir.
   """
    return lambda fname: base_dir / fname


def check_exists(fname: Path, file_description="File"):
    if not fname.exists():
        error_message = f"{file_description} required but not found: {fname}"
        log.error(error_message)
        exit(1)


class ProjectPaths:
    all_vars = {}

    def items(self):
        all_vars = vars(self)
        ProjectPaths.all_vars = {
            var: entry for var, entry in all_vars.items() if isinstance(entry, Path)
        }
        return ProjectPaths.all_vars.items()

    def dict(self):
        all_vars = vars(self)
        returned_vars = {}
        for var, entry in all_vars.items():
            if not isinstance(entry, Path):
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
        self.path_maker = path_fact(self.gram_dir)

        self.prg = self.path_maker("prg")
        self.ref = self.path_maker("original_reference.fasta")
        self.built_vcf = self.path_maker("build.vcf")
        self.report = self.path_maker("build_report.json")
        self.fm_index = self.path_maker("fm_index")
        self.cov_graph = self.path_maker("cov_graph")
        self.made_gram_dir = False

    def make_root(self):
        if not self.gram_dir.exists():
            log.debug("Creating gram directory:\n%s", self.gram_dir)
            self.gram_dir.mkdir()
            self.made_gram_dir = True

    def unmake_root(self):
        if self.gram_dir.exists() and self.made_gram_dir:
            self.gram_dir.rmdir()

    def _check_ref_and_make_ref_link(self, ref_path):
        resolved_ref_path = Path.resolve(Path(ref_path))
        check_exists(resolved_ref_path)
        if os.path.lexists(self.ref):
            self.ref.unlink()
        os.symlink(resolved_ref_path, self.ref)

    def _check_and_flatten_vcf_filenames(self, vcf: List[List[str]]):
        vcf_files = [Path(vcf_file) for arglist in vcf for vcf_file in arglist]

        for vcf_file in vcf_files:
            check_exists(vcf_file)
        self.input_vcfs = vcf_files

    def ready_ref_and_vcf(self, reference: str, vcfs: List[List[str]]):
        self._check_ref_and_make_ref_link(reference)
        self._check_and_flatten_vcf_filenames(vcfs)

    def raise_error(self, err_message):
        self.unmake_root()
        log.error(err_message)
        exit(1)


## Generates paths that will be shared between 'run' commands.
def _generate_run_paths(run_dir):

    run_dir = os.path.abspath(run_dir)

    # Quasimap paths
    run_path = path_fact(run_dir)

    paths = {
        "run_dir": run_dir,
        "build_dir": run_path(
            "build_dir"
        ),  #  This will be a symlink of the actual build dir
        "reads_dir": run_path(
            "reads"
        ),  # This will contain symlinks to the actual reads
        "quasimap_dir": run_path("quasimap_outputs"),
        "infer_dir": run_path("infer_outputs"),
        "discover_dir": run_path("discover_outputs"),
    }

    quasimap_path = path_fact(paths["quasimap_dir"])

    paths.update(
        {
            "allele_base_coverage": quasimap_path("allele_base_coverage.json"),
            "grouped_allele_counts_coverage": quasimap_path(
                "grouped_allele_counts_coverage.json"
            ),
            "allele_sum_coverage": quasimap_path("allele_sum_coverage"),
            "read_stats": quasimap_path("read_stats.json"),
        }
    )

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


## Make quasimap-related file and directory paths.
def generate_quasimap_paths(args):
    all_paths = _generate_project_paths(args.gram_dir)
    all_paths.update(_generate_run_paths(args.run_dir))

    if not os.path.exists(all_paths["run_dir"]):
        os.mkdir(all_paths["run_dir"])

    if not os.path.exists(all_paths["quasimap_dir"]):
        os.mkdir(all_paths["quasimap_dir"])

    # Make a reference to the gram_dir made in build. This avoids downstream commands
    # to require a --gram_dir argument.
    if os.path.lexists(all_paths["build_dir"]):
        os.unlink(all_paths["build_dir"])
    os.symlink(all_paths["gram_dir"], all_paths["build_dir"], target_is_directory=True)

    if os.path.exists(all_paths["reads_dir"]):
        log.warning(
            "Erasing reads directory {} and its contents to avoid accumulation".format(
                all_paths["reads_dir"]
            )
        )
        shutil.rmtree(all_paths["reads_dir"])

    os.mkdir(all_paths["reads_dir"])

    #  Make symlinks to the read files
    reads_files = [
        read_file for arglist in args.reads for read_file in arglist
    ]  # Flattens out list of lists.
    for read_file in reads_files:
        read_file = os.path.abspath(
            read_file
        )  # So that symlinking does not need absolute fpaths passed at command line.

        base = os.path.basename(read_file)
        target = os.path.join(all_paths["reads_dir"], base)
        if not os.path.exists(target):
            os.symlink(read_file, target)

    all_paths.update(
        {
            "report": os.path.join(all_paths["quasimap_dir"], "quasimap_report.json"),
            "reads_files": reads_files,
        }
    )

    return all_paths


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
