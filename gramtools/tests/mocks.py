class _MockVcfRecord:
    def __init__(self, pos, ref, alts, samples=[], chrom="JAC", filters=dict()):
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self.chrom = chrom

        if len(samples) == 0:
            self.samples = [
                {"GT": [1, 1]}
            ]  # Default to recording a ALT call for a single sample.
        else:
            self.samples = samples

        if len(filters) == 0:
            self.filters = {"PASS": ""}
        else:
            self.filters = filters

    def __repr__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        return (
            self.pos == other.pos
            and self.ref == other.ref
            and self.alts == other.alts
            and self.chrom == other.chrom
        )
