from pybedtools import BedTool, Interval
from enum import Enum, auto
from typing import List, Optional
import re
from pathlib import Path

from gramtools.commands.common import load_fasta

Intervals = List[Interval]

msa_like = re.compile("(msa|fa|fasta)$")


class BuildType(Enum):
    PRG = auto()
    MSA = auto()
    INVARIANT = auto()


def check_all_fnames_exist(intervals: Intervals) -> Optional[str]:
    for interval in intervals:
        fname = interval.name
        if not Path(fname).exists():
            return fname
        return None


class IntervalCollection:
    def __init__(self, bed_fname: str, fasta_fname: str, out_dirname: str):
        self.builders: List[IntervalBuilder] = list()
        input_collection = BedTool(bed_fname)
        missing_fname = check_all_fnames_exist(input_collection)
        if missing_fname is not None:
            raise ValueError(
                f"Error: {missing_fname} not found (specified in {bed_fname})"
            )
        for interval in input_collection:
            if msa_like.match(interval.name) is not None:
                build_type = BuildType.MSA
            else:
                build_type = BuildType.PRG
            out_fname = f"{out_dirname}/{Path(interval.name).stem}.bin"
            self.builders.append(IntervalBuilder(interval, build_type, out_fname))

        self.genome = load_fasta(fasta_fname)
        self.chrom_sizes = {name: len(seq) for name, seq in self.genome.items()}

        # TODO: add chrom dict to complement() and sort()
        invar_intervals = BedTool(self.intervals).complement()
        for i, invar_interval in enumerate(invar_intervals):
            out_fname = f"{out_dirname}/invariant_{i+1}.bin"
            self.builders.append(
                IntervalBuilder(invar_interval, BuildType.INVARIANT, out_fname)
            )
        print(self.intervals)

    @property
    def intervals(self):
        return [builder.interval for builder in self.builders]

    def dispatch_building(self):
        # Call build() on each interval builder
        # Can be multiprocessed
        pass

    def get_built_bed(self, out_fname):
        # Check each member interval has type PRG and an existing out_fname
        # Make set of intervals with added out_fname in col 4
        # Return bed. That bed can be written but also then used by concat_prg
        pass


class IntervalBuilder:
    def __init__(self, interval: Interval, build_type: BuildType, out_fname: str):
        self.chrom = interval.chrom
        self.start = interval.start
        self.end = interval.end
        self.out_fname = out_fname
        self.build_type = build_type
        self._in_fname = interval.name

    def build(self):
        # Run make_prg on MSA type
        # Copy existing prg
        pass

    @property
    def interval(self):
        return Interval(self.chrom, self.start, self.end, self.out_fname)


if __name__ == "__main__":
    import sys

    bed_fname = sys.argv[1]
    fasta_fname = ""
    out_dirname = "."
    IntervalCollection(bed_fname, fasta_fname, out_dirname)
