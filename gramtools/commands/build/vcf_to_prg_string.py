"""
Assumption: the vcf records are sorted by i)CHROM and then ii)position.

Behaviour:
    * (logged) For a given CHROM, records with overlapping or non-increasing POS are dropped.
    * (logged) Records without 'PASS' in FILTER are skipped.
    * Reference records with no variation are appended at the end of the prg string, in the order they appear
     in the ref file.
"""
from collections import OrderedDict, defaultdict
import logging
import sys

from pysam import VariantFile, VariantRecord
from Bio import SeqIO

sys.tracebacklimit = 0
logger = logging.getLogger("vcf_to_prg_string")
logger.setLevel(logging.WARNING)


def load_fasta(reference_file):
    ref_records = OrderedDict()
    # SeqIO: the id of the record is everything between ">" and the first space.
    for seq_record in SeqIO.parse(reference_file, "fasta"):
        ref_records[seq_record.id] = str(seq_record.seq)
    return ref_records


def nucleotide_to_integer(nucleotide):
    if nucleotide == "A":
        return 1
    elif nucleotide == "C":
        return 2
    elif nucleotide == "G":
        return 3
    elif nucleotide == "T":
        return 4
    else:
        raise ValueError("Did not receive a nucleotide ({A,C,G,T}")


def integer_to_nucleotide(integer):
    if integer == 1:
        return "A"
    elif integer == 2:
        return "C"
    elif integer == 3:
        return "G"
    elif integer == 4:
        return "T"
    else:
        return str(integer)


class ReferenceError(Exception):
    pass


class _Template_Vcf_to_prg(object):
    acceptable_modes = {"legacy", "normal"}
    NUM_BYTES = 4  # Number of bytes for each serialised integer

    def __init__(self):
        self.prg_vector: Map[str, List[int]] = defaultdict(list)
        self.num_sites = 0
        self.processed_refs = list()
        self.f_out = None
        self.skipped_records: int = 0

    def _check_record_ref(self, rec: VariantRecord):
        # Make sure the vcf record refers to a valid ref sequence ID
        if rec.chrom not in self.ref_records:
            raise ReferenceError(
                f"ref ID {rec.chrom} not found in reference file {self.ref_in}"
            )
        else:
            pos, length = rec.pos - 1, len(rec.ref)
            if self.ref_records[rec.chrom][pos : pos + length] != rec.ref:
                raise ReferenceError(
                    f"Vcf record REF sequence does not match ref ID {rec.chrom} "
                )

    def _record_to_prg_rep(self, site_marker, mode):
        """
        Go from a vcf record to its prg string representation:
            * (legacy) site markers flank the variant site
            * (normal) site marker at start, allele marker everywhere else.
        """
        prg_rep = [site_marker]
        allele_marker = site_marker + 1
        alleles = [[nucleotide_to_integer(nt) for nt in self.vcf_record.ref]]
        for alt_allele in self.vcf_record.alts:
            alleles.append([nucleotide_to_integer(nt) for nt in str(alt_allele)])
        if mode == "normal":
            for allele in alleles:
                prg_rep.extend(allele + [allele_marker])

        elif mode == "legacy":
            for i, allele in enumerate(alleles):
                if i < len(alleles) - 1:
                    prg_rep.extend(allele + [allele_marker])
                else:
                    prg_rep.extend(allele + [site_marker])
        return prg_rep

    def _get_ref_slice(self, chrom, start, end=0):
        ref_seq = self.ref_records[chrom]
        res = ref_seq[start - 1 :].upper()
        if end != 0:
            res = res[: end - start]
        res = list(map(nucleotide_to_integer, res))  # Convert to integers
        return res

    def _get_ref_records_with_no_variants(self):
        used_refs = set(self.processed_refs)
        diff = [rec_id for rec_id in self.ref_records if rec_id not in used_refs]
        if len(diff) > 0:
            for id in diff:
                seq = [nucleotide_to_integer(nt) for nt in self.ref_records[id]]
                yield id, seq


class Vcf_to_prg(_Template_Vcf_to_prg):
    """
    Args:
        mode : if 'legacy' a T/G SNP is built as '5T6G5',
               if 'normal' is built as '5T6G6'
    """

    def __init__(self, vcf_file, reference_file, prg_output_file, mode="normal"):
        super().__init__()
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
                    self.prg_vector[ref_chrom] += self._get_ref_slice(
                        ref_chrom, ref_pos
                    )
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
                self.prg_vector[ref_chrom] += self._get_ref_slice(
                    vcf_chrom, ref_pos, vcf_pos
                )
                ref_pos = vcf_pos

            self.prg_vector[ref_chrom] += self._record_to_prg_rep(cur_site_marker, mode)
            ref_pos += len(self.vcf_record.ref)
            cur_site_marker += 2

        # End condition: might have some invariant reference left
        if ref_chrom is not None:
            self.prg_vector[ref_chrom] += self._get_ref_slice(ref_chrom, ref_pos)
            self.processed_refs.append(vcf_chrom)

        # End condition #2 : might have records with no variants
        for ref_chrom, ref_seq in self._get_ref_records_with_no_variants():
            self.prg_vector[ref_chrom] += ref_seq
            self.processed_refs.append(ref_chrom)

    def _get_ints(self):
        contiguous = [self.prg_vector[ref_chrom] for ref_chrom in self.ref_records]
        return "".join(contiguous)

    def _get_string(self):
        prg_string = ""
        for ref_chrom in self.ref_records:
            for x in self.prg_vector[ref_chrom]:
                prg_string += integer_to_nucleotide(x)
        return prg_string

    def _write_bytes(self):
        """
        Write an integer
        By default, writes an **unsigned integer**
        Endianness: big or little
        :return:
        """
        with open(f"{self.f_out_prefix}", "wb") as f_out:
            for ref_chrom in self.ref_records:
                for integer in self.prg_vector[ref_chrom]:
                    f_out.write(integer.to_bytes(self.NUM_BYTES, "little"))

    def _write_string(self):
        # Construct it
        prg_string = ""
        for ref_chrom in self.ref_records:
            for x in self.prg_vector[ref_chrom]:
                prg_string += integer_to_nucleotide(x)
        with open(f"{self.f_out_prefix}.prg", "w") as f_out:
            f_out.write(prg_string)

    def _write_coordinates(self):
        with open(f"{self.f_out_prefix}_coords.tsv", "w") as genome_file:
            for ref_chrom, ref_seq in self.ref_records.items():
                line = f"{ref_chrom}\t{len(ref_seq)}\n"
                genome_file.write(line)


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
    converter._write_coordinates()
    logger.info(f"num variant sites in prg: {converter.num_sites}")
