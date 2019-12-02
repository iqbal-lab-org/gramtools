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

import collections
import vcf
from Bio import SeqIO

class _Template_Vcf_to_prg(object):
    """
    Defines utility functions for its descendants.
    Not for users.
    """
    acceptable_modes = {"legacy", "normal"}
    NUM_BYTES = 4 # Number of bytes for each serialised integer

    def __init__(self):
        self.prg_vector = []
        self.ref_records = collections.OrderedDict()
        self.num_sites = 0
        self.processed_ref = set()  # A set of ref records with variation
        self.f_out = None
        self.mode = None


    def _parse_reference(self, reference_file):
        """
        Read in the records in the fasta reference
        """
        # A note on SeqIO: the id of the record is everything between ">" and the first space.
        for seq_record in SeqIO.parse(reference_file, "fasta"):
            self.ref_records[seq_record.id] = seq_record.seq

    def _check_record_ref(self, chrom):
        """
        Make sure the vcf record refers to a valid ref sequence ID
        """
        ref_seq = self.ref_records.get(chrom, None)
        if ref_seq is None:
            raise ValueError(f"The ref ID {chrom} "
                             "in vcf file {self.vcf_in} was not found in reference file {self.ref_in}")

    def _record_to_prg_rep(self, vcf_record, site_marker):
        """
        Go from a vcf record to its prg string representation:
            * (legacy) site markers flank the variant site
            * (normal) site marker at start, allele marker everywhere else.
        """
        prg_rep = [site_marker]
        allele_marker = site_marker + 1
        alleles = [[integer_encode_nucleotide(nt) for nt in vcf_record.REF]]
        for alt_allele in vcf_record.ALT:
            alleles.append([integer_encode_nucleotide(nt) for nt in str(alt_allele)])
        if self.mode == "normal":
            for allele in alleles:
                prg_rep.extend(allele + [allele_marker])

        elif self.mode == "legacy":
            for i,allele in enumerate(alleles):
                if i < len(alleles) - 1:
                    prg_rep.extend(allele + [allele_marker])
                else:
                    prg_rep.extend(allele + [site_marker])
        return prg_rep




    def _get_invariant_portion(self, chrom, processed_pos, cur_pos):
        """
        Get the ref sequence in between records
        Assumes that `chrom` is valid, ie exists as a ref record.
        """
        ref_seq = self.ref_records[chrom]
        invariant = str(ref_seq[processed_pos - 1: cur_pos - 1]).upper()
        invariant = list(map(integer_encode_nucleotide, invariant)) # Convert to integers
        return invariant


    def _get_final_invariant(self, chrom, processed_pos):
        invariant = []
        ref_seq = self.ref_records[chrom]
        # The only case where we don't want to enter this is when the processed_pos has gone past the end of the ref record
        if len(ref_seq) >= processed_pos:
            invariant = ref_seq[processed_pos - 1:]

        invariant = list(map(integer_encode_nucleotide, invariant.upper())) # Convert to integers
        return invariant

    def _maybe_make_chrom_end(self, prev_chrom, cur_chrom, processed_pos):
        """
        When we change chromosome, we may need to add its invariant end to the prg_string.
        :return: Updated, or unchanged, chromo and current position in chromo.
        """
        if cur_chrom != prev_chrom:
            # Case: deal with end of chromosome
            if prev_chrom is not None:
                self.prg_vector.extend(self._get_final_invariant(prev_chrom, processed_pos))
                processed_pos = 1
            prev_chrom = cur_chrom
            self.processed_ref.add(cur_chrom)

        return prev_chrom, processed_pos


    def _maybe_add_preceding_invariant_sequence(self, chrom, cur_pos, processed_pos):
        """
        When we go to a new record, we may need to add invariant sequence preceding it to the prg_string.
        :return: Updated, or unchanged current position in chromo
        """
        if cur_pos > processed_pos:
            #  Case: add invariant reference in between variant records
            self.prg_vector += self._get_invariant_portion(chrom, processed_pos, cur_pos)


    def _get_ref_records_with_no_variants(self):
        """
        Writes each record with no variants to the prg file
        """
        ref_no_recs = []
        all_records = set(self.ref_records.keys())
        diff = all_records.difference(self.processed_ref)
        diff = [r for r in self.ref_records.keys() if r in diff] # Conserve the initial ordering.
        if len(diff) > 0:
            for id in diff:
                seq = [integer_encode_nucleotide(nt) for nt in self.ref_records[id]]
                ref_no_recs.append(seq)
        return ref_no_recs



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
        self.mode = mode

        self.vcf_in = vcf_file
        self.ref_in = reference_file
        self.f_out_prefix = prg_output_file
        self._parse_reference(reference_file) # Populates the ref_records
        self._make_prg();


    def _make_prg(self):
        """
        Constructs the prg string.
        """
        cur_site_marker = 5
        processed_pos = 1  # 1-based
        prev_chrom = None

        #  Parse vcf
        f_in = open(self.vcf_in)
        vcf_in = vcf.Reader(f_in)

        try:
            record = next(vcf_in, None)
        except Exception as e: # PyVCF raises this
            print(f'ERROR: Could not read vcf file {self.vcf_in} due to: "{e}".'
                  'Is the VCF properly formatted?')
            raise

        while record != None:
            cur_pos = record.POS

            self._check_record_ref(record.CHROM)

            prev_chrom, processed_pos = self._maybe_make_chrom_end(prev_chrom, record.CHROM, processed_pos)

            # Skip records which are not sorted, or all others present at same position
            if cur_pos < processed_pos:
                record = next(vcf_in, None)
                continue

            self._maybe_add_preceding_invariant_sequence(record.CHROM, cur_pos, processed_pos)

            # Convert vcf record to its prg representation
            self.prg_vector.extend(self._record_to_prg_rep(record, cur_site_marker))
            processed_pos = cur_pos + len(record.REF)

            cur_site_marker += 2
            self.num_sites += 1

            record = next(vcf_in, None)

        # End condition: might have some invariant reference left
        if prev_chrom is not None:
            self.prg_vector.extend(self._get_final_invariant(prev_chrom, processed_pos))

        # End condition #2 : might have records with no variants
        for rec in self._get_ref_records_with_no_variants():
            self.prg_vector.extend(rec)

        f_in.close()

    def _write_bytes(self):
        """
        Write an integer
        By default, writes an **unsigned integer**
        Endianness: big or little
        :return:
        """
        with open(f'{self.f_out_prefix}', "wb") as f_out:
            for integer in self.prg_vector:
                f_out.write(integer.to_bytes(self.NUM_BYTES, "little"))

    def _write_string(self):
        # Construct it
        prg_string = ""
        for x in self.prg_vector:
            prg_string += nucleotide_decode_integer(x)
        with open(f'{self.f_out_prefix}.prg', "w") as f_out:
            f_out.write(str(prg_string));


# May run standalone
if __name__ == "__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('vcf', type=str, help='vcf file in')
    parser.add_argument('ref', type=str, help='reference file in')
    parser.add_argument('--outfile', type=str, help='path to prg file out (default: file "prg" in cwd)',
                        default=os.path.join(os.path.abspath(os.path.curdir), "prg"))
    parser.add_argument('--mode', type=str, help='Run in legacy or normal mode',
                        choices=["legacy","normal"],
                        default="normal")

    args = parser.parse_args()

    converter = Vcf_to_prg(args.vcf, args.ref, args.outfile, args.mode)
    converter._write_bytes()
    converter._write_string()
    print(f'num variant sites in prg: {converter.num_sites}')
