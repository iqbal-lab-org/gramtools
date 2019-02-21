class GenomeRegions:
    """Variantblocks container class.

    Contains entire genome split into variant blocks.
    """

    def __init__(self):
        self._regions = []
        self._region_idx = {}  # Maps a `_GenomeRegion` to its index in `_regions`

    def __iter__(self):
        for region in self._regions:
            yield region

    def __getitem__(self, idx):
        return self._regions[idx]

    def range(self, start_region, reverse=False):
        start_region_idx = self._region_idx[start_region]
        if not reverse:
            idx_range = range(start_region_idx + 1, len(self._regions))
        else:
            idx_range = range(start_region_idx - 1, -1, -1)

        return (self._regions[idx] for idx in idx_range)

    def add_region(self, alleles, variant_site_marker=None):
        """Add a region of the genome."""
        region = _GenomeRegion(alleles, variant_site_marker)
        self._regions.append(region)
        self._region_idx[region] = len(self._regions) - 1

    def __str__(self):
        return str(''.join(str(region) for region in self._regions))

    def __repr__(self):
        return str(self)


class _GenomeRegion:
    """A region of the genome: either a variant site or a non-variant site."""

    def __init__(self, str_alleles, variant_site_marker=None):
        self.alleles = [tuple(a) for a in str_alleles]
        self.min_allele_len = min(len(i) for i in self.alleles)
        self.max_allele_len = max(len(i) for i in self.alleles)
        self.variant_site_marker = variant_site_marker

    @property
    def is_variant_site(self):
        return self.variant_site_marker is not None

    def __str__(self):
        return '[{alleles}]'.format(alleles='|'.join(str(x) for x in self.alleles))

    def __repr__(self):
        return str(self)
