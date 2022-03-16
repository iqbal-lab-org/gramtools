"""
Assumption: the vcf records are sorted by i)CHROM and then ii)position.

Behaviour:
    * (logged) For a given CHROM, records with overlapping or non-increasing POS are dropped.
    * (logged) Records without 'PASS' in FILTER are skipped.
    * Reference records with no variation are appended at the end of the prg string, in the order they appear
     in the ref file.
"""
from collections import defaultdict
import logging
import sys
from typing import List, Tuple

from pysam import VariantFile, VariantRecord

from gramtools import ENDIANNESS, BYTES_PER_INT
from gramtools.commands.common import load_fasta, int_to_bytes, integer_to_nucleotide

logger = logging.getLogger("vcf_to_prg_string")
logger.setLevel(logging.WARNING)


class ReferenceError(Exception):
    pass


class Vcf_to_prg(object):
    """
    Args:
        mode : if 'legacy' a T/G SNP is built as '5T6G5',
               if 'normal' is built as '5T6G6'
    """

    acceptable_modes = {"legacy", "normal"}

    def __init__(self, vcf_file, reference_file, prg_output_file, mode="normal"):
        self.prg_bytes: Map[str, bytearray] = defaultdict(bytearray)
        self.num_sites = 0
        self.processed_refs = list()
        self.skipped_records: int = 0

        if mode not in self.acceptable_modes:
            raise ValueError(f"Mode not in {self.acceptable_modes}")

        self.f_out_prefix = prg_output_file
        self.vcf_in = VariantFile(vcf_file).fetch()
        self.ref_in = reference_file

        self.ref_records = load_fasta(reference_file)

        self._make_prg(mode)
        if self.skipped_records > 0:
            logger.warning(
                f"Skipped {self.skipped_records} because of no 'PASS' in their FORMAT column"
            )

    def _from_bytes(self, to_convert: bytes) -> int:
        return int.from_bytes(to_convert, ENDIANNESS)

    def _get_ref_slice(self, chrom, start, end=0) -> bytearray:
        if end == 0:
            res = self.ref_records[chrom][start - 1 :]
        else:
            res = self.ref_records[chrom][start - 1 : end - 1]
        return bytearray(b"".join(map(int_to_bytes, res)))

    def _check_record_ref(self, rec: VariantRecord) -> None:
        # Make sure the vcf record refers to a valid ref sequence ID
        if rec.chrom not in self.ref_records:
            raise ReferenceError(
                f"ref ID {rec.chrom} not found in reference file {self.ref_in}"
            )
        else:
            pos, length = rec.pos - 1, len(rec.ref)
            if self.ref_records[rec.chrom][pos : pos + length].upper() != rec.ref:
                raise ReferenceError(
                    f"Vcf record REF sequence does not match ref ID {rec.chrom} "
                )

    def _record_to_prg_rep(self, site_marker: int, mode: str) -> bytearray:
        """
        Go from a vcf record to its prg string representation:
            * (legacy) site markers flank the variant site
            * (normal) site marker at start, allele marker everywhere else.
        """
        # Add initial site marker + ref allele
        prg_rep: bytearray = bytearray(int_to_bytes(site_marker))
        prg_rep += b"".join(map(int_to_bytes, self.vcf_record.ref))

        # Add allele markers + alt alleles
        allele_marker = site_marker + 1
        prg_rep += int_to_bytes(allele_marker)
        num_alts = len(self.vcf_record.alts)
        for i, alt_allele in enumerate(self.vcf_record.alts):
            prg_rep += b"".join(map(int_to_bytes, str(alt_allele)))
            pushed_marker = allele_marker
            if mode == "legacy" and i == num_alts - 1:
                pushed_marker -= 1
            prg_rep += int_to_bytes(pushed_marker)
        return prg_rep

    def _get_ref_records_with_no_variants(self) -> Tuple[str, bytearray]:
        diff = set(self.ref_records).difference(set(self.processed_refs))
        for id in diff:
            seq = bytearray(b"".join(map(int_to_bytes, self.ref_records[id])))
            yield id, seq

    def next_rec(self):
        try:
            while True:
                self.vcf_record = next(self.vcf_in)
                if "PASS" not in self.vcf_record.filter:
                    self.skipped_records += 1
                    continue
                return True
        except StopIteration:
            return False
        except Exception as e:
            logger.error(
                f'ERROR: Could not read vcf file {self.vcf_in} due to: "{e}".'
                "Is the VCF properly formatted?"
            )
            raise

    def _make_prg(self, mode):
        ref_chrom = None
        cur_site_marker = 5

        while self.next_rec():
            self.num_sites += 1
            vcf_pos = self.vcf_record.pos
            vcf_chrom = self.vcf_record.chrom
            self._check_record_ref(self.vcf_record)

            if vcf_chrom != ref_chrom:
                if ref_chrom is not None:
                    self.prg_bytes[ref_chrom] += self._get_ref_slice(ref_chrom, ref_pos)
                    self.processed_refs.append(ref_chrom)
                ref_pos = 1
                ref_chrom = vcf_chrom

            # Skip records which are not sorted, or all others present at same position
            if vcf_pos < ref_pos:
                logger.warning(
                    f"Skipped record at pos {vcf_pos}, because previous record led to pos {ref_pos}"
                )
                continue

            if vcf_pos > ref_pos:
                self.prg_bytes[ref_chrom] += self._get_ref_slice(
                    vcf_chrom, ref_pos, vcf_pos
                )
                ref_pos = vcf_pos

            self.prg_bytes[ref_chrom] += self._record_to_prg_rep(cur_site_marker, mode)
            ref_pos += len(self.vcf_record.ref)
            cur_site_marker += 2

        # End condition: might have some invariant reference left
        if ref_chrom is not None:
            self.prg_bytes[ref_chrom] += self._get_ref_slice(ref_chrom, ref_pos)
            self.processed_refs.append(vcf_chrom)

        # End condition #2 : might have records with no variants
        for ref_chrom, ref_seq in self._get_ref_records_with_no_variants():
            self.prg_bytes[ref_chrom] += ref_seq
            self.processed_refs.append(ref_chrom)

    def _get_ints(self) -> List[int]:
        contiguous_ints = list()
        for ref_chrom in self.ref_records:
            chrom_bytes = self.prg_bytes[ref_chrom]
            for inset in range(0, len(chrom_bytes), BYTES_PER_INT):
                contiguous_ints.append(
                    self._from_bytes(chrom_bytes[inset : inset + BYTES_PER_INT])
                )
        return contiguous_ints

    def _get_string(self):
        prg_ints = self._get_ints()
        return "".join(map(integer_to_nucleotide, prg_ints))

    def _write_bytes(self):
        bytes = bytearray()
        for ref_chrom in self.ref_records:
            bytes += self.prg_bytes[ref_chrom]
        with open(f"{self.f_out_prefix}", "wb") as f_out:
            f_out.write(bytes)

    def _write_string(self):
        prg_string = self._get_string()
        with open(f"{self.f_out_prefix}.prg", "w") as f_out:
            f_out.write(prg_string)


# May run standalone
if __name__ == "__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("vcf", type=str, help="vcf file in")
    parser.add_argument("ref", type=str, help="reference file in")
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help='path to prg file out (default: file "prg" in cwd)',
        default=os.path.join(os.path.abspath(os.path.curdir), "prg"),
    )
    parser.add_argument(
        "--mode",
        type=str,
        help="Run in legacy or normal mode",
        choices=["legacy", "normal"],
        default="normal",
    )

    args = parser.parse_args()

    converter = Vcf_to_prg(args.vcf, args.ref, args.outfile, args.mode)
    converter._write_bytes()
    converter._write_string()
    logger.info(f"num variant sites in prg: {converter.num_sites}")
