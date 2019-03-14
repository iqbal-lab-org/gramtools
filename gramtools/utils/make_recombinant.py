## @file Make new fasta (recombinant) from fasta + vcf.

from .. import prg_local_parser
import vcf
import random
import os

DNA = {'A', 'C', 'G', 'T'}
REF_READ_CHUNK_SIZE = 1000000

## An interface to the module for making new fasta (recombinant) from fasta + vcf.
class Recombinator(object):
    def __init__(self, reference_fpath, vcf_fpath, output_fpath, random_selection = True, seed = None, offset = 1):
        self.reference_fpath = reference_fpath
        self.vcf_fpath = vcf_fpath
        self.output_fpath = output_fpath
        self.random_selection = random_selection
        self.seed = seed
        # This is used to signal that we do not start counting from 0 when traversing the original fasta.
        # We use 1 as a default value to go from 1-based VCF coords to 0-based `enumerate()`
        self.offset = offset

    # 'static' variable
    picked_alleles = []

    def run(self):
        fhandle = open(self.vcf_fpath)
        vcf_reader = vcf.Reader(fhandle)
        description = "Recombinant built from ref and vcf"
        _make_recombinant(self.reference_fpath, vcf_reader, self.output_fpath, description, self.random_selection, self.seed, self.offset)
        if len(Recombinator.picked_alleles) > 0:
            picked_file = os.path.realpath(os.path.dirname(self.output_fpath))
            picked_file = os.path.join(picked_file, "picked_alleles")
            with open(picked_file, "w") as f:
                f.write("\t".join(Recombinator.picked_alleles))
        fhandle.close()

## Compute a new reference from an base reference and a set of vcf records.
# @param reference the 'base' reference: non-variant sites of the prg + first allele of each variant site.
# @param vcf_reader set of vcf records describing variants inferred against the old reference.
# Default random_selection to False because can be used as is for 'genotyping' mode.
def _make_recombinant(reference_fpath, vcf_reader, output_fpath, description, random_selection = False, seed = None, offset = 1):

    if random_selection and seed is not None:
        random.seed(seed)

    inferred_reference_length = 0

    reference = _refParser(reference_fpath)
    recombinant = prg_local_parser.FastaWriter(output_fpath, description)

    record = next(vcf_reader, None)

    start_ref_idx = -1
    end_ref_idx = -1

    for idx, x in enumerate(reference):
        within_replaced_ref = start_ref_idx < idx <= end_ref_idx
        if within_replaced_ref:
            continue

        vcf_finished = record is None
        if vcf_finished:
            recombinant.append(x)
            inferred_reference_length += 1
            continue

        # This update only actually occurs once we have passed a record
        start_ref_idx = record.POS - offset
        end_ref_idx = start_ref_idx + len(record.REF) - 1

        at_allele_insert_idx = idx == start_ref_idx

        # Let's insert the genotyped allele
        if at_allele_insert_idx:
            sample = record.samples[0]
            genotype = sample.gt_alleles

            if random_selection is False: # Genotyping mode
                # Case: no genotype. Let's pick REF allele.
                if set(genotype) == {None}:
                    allele_index = 0

                # Case: single haploid genotype. Pick that allele.
                # OR Case: there are >1 calls. Just pick the first.
                else:
                        allele_index = int(genotype[0])

            else: # Random selection mode
                allele_index = random.randrange(len(record.ALT) + 1) # Add 1 to choose between REF and all ALTs
                Recombinator.picked_alleles.append(str(allele_index))

            if allele_index == 0:
                template = str(record.REF)
            else:
                template = str(record.ALT[allele_index - 1])

            for alt_char in template:
                recombinant.append(alt_char)

            inferred_reference_length += len(template)

            record = next(vcf_reader, None)
            continue

        # Add the non-variant regions
        recombinant.append(x)
        inferred_reference_length += 1

    recombinant.close()

    return inferred_reference_length


## Generator to parse fasta ref sequence in blocks.
def _refParser(ref_fpath, chunk_size = REF_READ_CHUNK_SIZE):
    with open(ref_fpath) as fhandle:
        next(fhandle)  # Skip header line
        while True:
            chars = fhandle.read(chunk_size)
            if len(chars) == 0:
                return

            for char in chars:
                if char.upper() not in DNA:
                    continue
                yield char


