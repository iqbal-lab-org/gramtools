from enum import Enum, auto
from typing import List, Dict, Optional
import re
import shutil
from pathlib import Path
import logging

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
from pybedtools import BedTool, Interval
from make_prg.from_msa.prg_builder import PrgBuilder
from make_prg.prg_encoder import PrgEncoder, ENDIANNESS, BYTES_PER_INT

from gramtools.commands.common import load_fasta, Chroms
from gramtools.commands.build.command_setup import MSA_EXTS
from gramtools.commands import report

Intervals = List[Interval]

msa_like = re.compile(MSA_EXTS)

log = logging.getLogger("gramtools")


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
    """
    Manages a list of `IntervalBuilder`s
    """

    def __init__(
        self, bed_fname: str, fasta_fname: str, coords_fname: str, out_dirname: str
    ):
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

        self.chrom_seqs: Chroms = load_fasta(fasta_fname)
        self.coords_fname = coords_fname

        invar_intervals = BedTool(self.intervals).complement(g=self.coords_fname)
        for i, invar_interval in enumerate(invar_intervals):
            out_fname = f"{out_dirname}/invariant_{i+1}.bin"
            new_builder = IntervalBuilder(
                invar_interval, BuildType.INVARIANT, out_fname
            )
            new_builder.add_seq(self.chrom_seqs)
            self.builders.append(new_builder)
        self.chrom_seqs = {}

    @property
    def intervals(self):
        return [builder.interval for builder in self.builders]

    def build(self):
        for builder in self.builders:
            builder.build()

    def get_built_bed(self):
        new_intervals = self.intervals
        check_all_fnames_exist(new_intervals)
        result = BedTool(new_intervals)
        result = result.sort(g=self.coords_fname)
        return result


class IntervalBuilder:
    """
    Responsible for building one prg
    """

    def __init__(self, interval: Interval, build_type: BuildType, out_fname: str):
        self.chrom = interval.chrom
        self.start = interval.start
        self.end = interval.end
        self.out_fname = out_fname
        self.build_type = build_type
        self._in_fname = interval.name

    def encode_and_write_prg(self, built_prg: PrgBuilder):
        prg_encoder = PrgEncoder()
        prg_ints = prg_encoder.encode(built_prg.prg)
        with open(self.out_fname, "wb") as fhandle_out:
            prg_encoder.write(prg_ints, fhandle_out)

    def add_seq(self, chrom_seqs: Chroms):
        self.sequence = chrom_seqs[self.chrom][self.start : self.end]

    def build(self):
        if self.build_type is BuildType.PRG:
            log.debug(f"Copying already-build prg {self._in_fname}")
            shutil.copy(self._in_fname, self.out_fname)
        elif self.build_type is BuildType.MSA:
            log.debug(f"Building variant prg from MSA {self._in_fname}")
            built_prg = PrgBuilder(self._in_fname)
            self.encode_and_write_prg(built_prg)
        else:
            log.debug(f"Building invariant prg for region {self.start}-{self.end}")
            records = [SeqRecord(Seq(self.sequence.upper()), id="invar_seq")]
            alignment = MultipleSeqAlignment(records)
            built_prg = PrgBuilder("_", alignment=alignment)
            self.encode_and_write_prg(built_prg)

    @property
    def interval(self):
        return Interval(self.chrom, self.start, self.end, self.out_fname)


class Record:
    def __init__(self, translation: int, count: int):
        self.translation = translation
        self.count = count


PRG_Ints = List[int]


class PRGAggregationError(Exception):
    pass


class PRGDecodeError(Exception):
    pass


class PRGAggregator:
    """
    Rescales variant site markers across multiple PRGs so that they are consistent
    E.g. two different PRGs will use '5' and '6' to delimit their first sites. 
    Aggregator makes sure that does not happen.
    """

    def __init__(self):
        self.translations: Dict[str, Dict[int, Record]] = dict()
        self.next_allocated = 5

    def translate(self, ID: str, marker: int) -> int:
        if ID not in self.translations:
            self.translations[ID] = dict()

        if marker <= 4:
            raise PRGAggregationError(f"Marker {marker} is not >4")

        local_table = self.translations[ID]
        if marker % 2 == 0:
            site_ID = marker - 1
            if site_ID not in local_table:
                raise PRGAggregationError(
                    f"Error: {marker}'s site number {marker - 1} has never been seen"
                )
            return local_table[site_ID].translation + 1

        else:
            if marker in local_table:
                record = local_table[marker]
                record.count += 1
                if record.count > 2:
                    raise PRGAggregationError(
                        f"Error: {marker} site number present >2 times in local PRG {ID}"
                    )
                else:
                    # Legacy format support: converts ending odd marker to even marker
                    return local_table[marker].translation + 1
            local_table[marker] = Record(self.next_allocated, 1)
            self.next_allocated += 2

            return local_table[marker].translation


def get_aggregated_prgs(agg: PRGAggregator, intervals: Intervals) -> PRG_Ints:
    rescaled_prg_ints: PRG_Ints = list()
    for interval in intervals:
        in_fname = interval.name
        prg_name = Path(in_fname).stem
        fhandle_in = open(in_fname, "rb")
        all_bytes = fhandle_in.read()
        fhandle_in.close()
        for pos in range(0, len(all_bytes), BYTES_PER_INT):
            int_bytes = all_bytes[pos : pos + BYTES_PER_INT]

            decoded_int = int.from_bytes(int_bytes, ENDIANNESS)
            if decoded_int <= 0:
                raise PRGDecodeError(f"PRG marker {decoded_int} should be > 0")
            elif decoded_int <= 4:
                rescaled_prg_ints.append(decoded_int)
            else:
                rescaled_int = agg.translate(prg_name, decoded_int)
                rescaled_prg_ints.append(rescaled_int)
    log.info(f"Total length of built prg: {len(rescaled_prg_ints)}")
    log.info(f"Total number of sites: {(agg.next_allocated - 3) // 2 - 1}\n")
    return rescaled_prg_ints


@report.with_report
def build_from_msas(report, action, build_paths, args):
    log.info("Building prg from prgs in {args.prgs_bed}")
    ic = IntervalCollection(
        args.prgs_bed,
        args.reference,
        build_paths.coords_file,
        build_paths.built_prg_dirname,
    )
    ic.build()
    built_intervals = ic.get_built_bed()
    built_intervals.saveas(build_paths.built_prg_bed)
    agg = PRGAggregator()
    rescaled_prg_ints: PRG_Ints = get_aggregated_prgs(agg, built_intervals)
    with open(build_paths.prg, "wb") as fhandle_out:
        PrgEncoder.write(rescaled_prg_ints, fhandle_out)
