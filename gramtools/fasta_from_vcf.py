## @file Make new fasta from fasta + vcf.

import vcf

DNA = {'A', 'C', 'G', 'T', 'N'}
REF_READ_CHUNK_SIZE = 1000000


def make_new_ref_using_vcf(original_fasta, vcf_file, out_fasta):
    """
    Assumption: the vcf records are sorted in same order as the original_fasta
    :param original_fasta:
    :param vcf_file:
    :param out_fasta:
    :return:
    """
    in_parse = _refParser(original_fasta)
    out_parse = _FastaWriter(out_fasta)

    vcf_Reader = vcf.Reader(open(vcf_file))

    vcf_record = next(vcf_Reader, None)
    start_pos = vcf_record.POS - 1
    end_pos = start_pos + len(vcf_record.REF) - 1

    inside_record = False
    start = True

    for fasta_char in in_parse:
        if fasta_char.header:
            reference_index = -1
            if not start:
                out_parse.commit('\n', is_used = fasta_char.use)
            else:
                start = False
        elif fasta_char.use:
            reference_index += 1
            if reference_index == start_pos: # Commit the chosen ALT
                inside_record = True

                # Pick the ALT based on GT column
                sample = vcf_record.samples[0]
                genotype = sample.gt_alleles

                # Case: no genotype. Let's pick REF allele.
                if set(genotype) == {None}:
                    chosen_pos = 0

                #Case: haploid, or >1 calls. Pick the first (is most likely haploid from infer).
                else:
                    chosen_pos = int(genotype[0])

                if chosen_pos == 0:
                    chosen_allele = vcf_record.REF
                else:
                    chosen_allele = vcf_record.ALT[chosen_pos - 1]

                for alt_char in str(chosen_allele):
                    out_parse.commit(alt_char)

            if reference_index == end_pos: # Get next record
                inside_record = False
                vcf_record = next(vcf_Reader, None)
                start_pos = vcf_record.POS - 1
                end_pos = start_pos + len(vcf_record.REF) - 1

            if inside_record:
                continue

        out_parse.commit(fasta_char.char, is_used = fasta_char.use)

    out_parse.close()


# Writes fasta characters in blocks
class _FastaWriter():
    max_load = REF_READ_CHUNK_SIZE
    line_length = 60

    def __init__(self, outfilename):
        self.load = []
        self.tally = 0
        self.fhandle = open(outfilename, "w")

    def commit(self, char, is_used = True):
        if is_used:
            char = char.upper()
        self.load.append(char)
        if is_used:
            self.tally += 1
        else:
            self.tally = 0

        if self.tally == _FastaWriter.line_length:
            self.load.append('\n')
            self.tally = 0
        if len(self.load) > _FastaWriter.max_load:
            self.fhandle.write("".join(self.load))
            self.load = []

    def close(self):
        if len(self.load) > 0:
            self.fhandle.write("".join(self.load))
        self.fhandle.close()

class Fasta_Char():
    valid = {"use", "header", "write"}

    def __init__(self, char, mode):
        self.char = char

        if mode not in Fasta_Char.valid:
            raise ValueError("Must initialise with mode in {}".format(str(Fasta_Char.valid)))

        if mode == "use":
            self.use = True
            self.header = False

        elif mode == "header":
            self.use = False
            self.header = True

        else:
            self.use = False
            self.header = False

## Generator to parse fasta ref sequence in blocks.
def _refParser(ref_fpath, chunk_size=REF_READ_CHUNK_SIZE):
    header = False
    with open(ref_fpath) as fhandle:
        while True:
            chars = fhandle.read(chunk_size)

            if len(chars) == 0:
                return

            for char in chars:
                if char == ">":
                    header = True
                    fasta_char = Fasta_Char(char, mode="header")

                elif char == '\n':
                    if header:
                        header = False
                        fasta_char = Fasta_Char(char, mode="write")
                    else: # Skip newlines if they are not part of the header
                        continue

                else:
                    if header:
                        fasta_char = Fasta_Char(char, mode="write")
                    else:
                        if char.upper() not in DNA:
                            raise ValueError("Found an invalid character: {}".format(char))
                        fasta_char = Fasta_Char(char, mode="use")

                yield fasta_char

