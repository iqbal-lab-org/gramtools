"""
Expects:
    A well formed vcf, here meaning no overlapping records- as output by vcf_clusterer module for eg.

Outputs:
    A prg string of the vcf, which `gramtools build` can use.

**Assumptions** (This module will NOT WORK AS YOU MAY WANT if NOT MET):
    * The vcf records are sorted by i)CHROM and then ii)position

Behaviour:
    * Reference records with no variation are all appended at the end of the prg string, in the order they first appear
     in the ref file.
    * (For a given CHROM,) Records with POS not strictly increasing are dropped.
"""
from collections import OrderedDict
import vcf
from Bio import SeqIO


def load_fasta(reference_file):
    ref_records = OrderedDict()
    # SeqIO: the id of the record is everything between ">" and the first space.
    for seq_record in SeqIO.parse(reference_file, "fasta"):
        ref_records[seq_record.id] = seq_record.seq
    return ref_records


def integer_encode_nucleotide(nucleotide):
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


def nucleotide_decode_integer(integer):
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


class _Template_Vcf_to_prg(object):
    acceptable_modes = {"legacy", "normal"}
    NUM_BYTES = 4  # Number of bytes for each serialised integer

    def __init__(self):
        self.prg_vector: List[int] = []
        self.num_sites = 0
        self.processed_refs = list()
        self.f_out = None

    def _check_record_ref(self, chrom):
        # Make sure the vcf record refers to a valid ref sequence ID
        if chrom not in self.ref_records:
            raise ValueError(
                f"The ref ID {chrom} "
                f"in vcf file {self.vcf_in} was not found in reference file {self.ref_in}"
            )

    def _record_to_prg_rep(self, site_marker, mode):
        """
        Go from a vcf record to its prg string representation:
            * (legacy) site markers flank the variant site
            * (normal) site marker at start, allele marker everywhere else.
        """
        prg_rep = [site_marker]
        allele_marker = site_marker + 1
        alleles = [[integer_encode_nucleotide(nt) for nt in self.vcf_record.REF]]
        for alt_allele in self.vcf_record.ALT:
            alleles.append([integer_encode_nucleotide(nt) for nt in str(alt_allele)])
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
        res = str(ref_seq[start - 1 :]).upper()
        if end != 0:
            res = res[: end - start]
        res = list(map(integer_encode_nucleotide, res))  # Convert to integers
        return res

    def _get_ref_records_with_no_variants(self):
        used_refs = set(self.processed_refs)
        diff = [rec_id for rec_id in self.ref_records.keys() if rec_id not in used_refs]
        if len(diff) > 0:
            for id in diff:
                seq = [integer_encode_nucleotide(nt) for nt in self.ref_records[id]]
                yield id, seq


class Vcf_to_prg(_Template_Vcf_to_prg):
    """
    The mode parameter, with default 'normal', defines how variant sites are built.
    If 'legacy', a T/G SNP is built as '5T6G5'
    If 'normal', a T/G SNP is built as '5T6G6'
    """

    def __init__(self, vcf_file, reference_file, prg_output_file, mode="normal"):
        super().__init__()
        if mode not in self.acceptable_modes:
            raise ValueError(f"Mode not in {self.acceptable_modes}")

        self.f_in = open(vcf_file)
        self.f_out_prefix = prg_output_file
        self.vcf_in = vcf.Reader(self.f_in)
        self.ref_in = reference_file

        self.ref_records = load_fasta(reference_file)
        self.vcf_record = self.next_rec()
        self._make_prg(mode)

    def next_rec(self):
        try:
            vcf_record = next(self.vcf_in, None)
        except Exception as e:  # PyVCF raises this
            print(
                f'ERROR: Could not read vcf file {self.vcf_in} due to: "{e}".'
                "Is the VCF properly formatted?"
            )
            self.f_in.close()
            raise
        return vcf_record

    def _make_prg(self, mode):
        ref_chrom = None
        cur_site_marker = 5

        while self.vcf_record != None:
            self.num_sites += 1
            vcf_pos = self.vcf_record.POS
            vcf_chrom = self.vcf_record.CHROM
            self._check_record_ref(vcf_chrom)

            if vcf_chrom != ref_chrom:
                if ref_chrom is not None:
                    self.prg_vector.extend(self._get_ref_slice(ref_chrom, ref_pos))
                    self.processed_refs.append(ref_chrom)
                ref_pos = 1
                ref_chrom = vcf_chrom

            # Skip records which are not sorted, or all others present at same position
            if vcf_pos < ref_pos:
                self.vcf_record = self.next_rec()
                continue

            if vcf_pos > ref_pos:
                self.prg_vector += self._get_ref_slice(vcf_chrom, ref_pos, vcf_pos)
                ref_pos = vcf_pos

            self.prg_vector.extend(self._record_to_prg_rep(cur_site_marker, mode))
            ref_pos += len(self.vcf_record.REF)
            cur_site_marker += 2

            self.vcf_record = self.next_rec()

        # End condition: might have some invariant reference left
        if ref_chrom is not None:
            self.prg_vector.extend(self._get_ref_slice(ref_chrom, ref_pos))
            self.processed_refs.append(vcf_chrom)

        # End condition #2 : might have records with no variants
        for ref_chrom, ref_seq in self._get_ref_records_with_no_variants():
            self.prg_vector.extend(ref_seq)
            self.processed_refs.append(ref_chrom)

        self.f_in.close()

    def _write_bytes(self):
        """
        Write an integer
        By default, writes an **unsigned integer**
        Endianness: big or little
        :return:
        """
        with open(f"{self.f_out_prefix}", "wb") as f_out:
            for integer in self.prg_vector:
                f_out.write(integer.to_bytes(self.NUM_BYTES, "little"))

    def _write_string(self):
        # Construct it
        prg_string = ""
        for x in self.prg_vector:
            prg_string += nucleotide_decode_integer(x)
        with open(f"{self.f_out_prefix}.prg", "w") as f_out:
            f_out.write(str(prg_string))


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
    print(f"num variant sites in prg: {converter.num_sites}")
