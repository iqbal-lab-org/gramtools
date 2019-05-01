import vcf

class _Sample:
    def __init__(self,genotype):
        self.gt_alleles = genotype


class _MockVcfRecord:
    def __init__(self, POS, REF, ALT, samples=[]):
        self.POS = POS
        self.REF = REF
        self.ALT = [vcf.model._Substitution(x) for x in ALT]

        if len(samples) == 0:
            self.samples = [_Sample(["1","1"])] # Default to recording a ALT call for a single sample.
        else:
            self.samples = samples

    def __repr__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        return self.POS == other.POS \
               and self.REF == other.REF \
               and self.ALT == other.ALT

