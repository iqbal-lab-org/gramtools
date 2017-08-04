class GenomeRegions:
    """Variantblocks container class.

    Contains entire genome split into variant blocks.
    """

    def __init__(self):
        self._regions = []
        self._region_idx = {}

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

        for idx in idx_range:
            region = self._regions[idx]
            yield region

    def add_region(self, region, variant_site_marker=None):
        """Add a region of the genome."""
        if not isinstance(region, _GenomeRegion):
            region = _GenomeRegion(region, variant_site_marker)

        self._regions.append(region)
        self._region_idx[region] = len(self._regions) - 1

    def __str__(self):
        return ''.join(str(region) for region in self._regions)


class _GenomeRegion:
    """A region of the genome, either a variant site or a non-variant site."""

    def __init__(self, alleles, variant_site_marker=None):
        # TODO: why is alleles copied here?
        self.alleles = list(alleles)
        self.min_allele_len = min(len(i) for i in alleles)
        self.max_allele_len = max(len(i) for i in alleles)
        self.variant_site_marker = variant_site_marker

    @property
    def is_variant_site(self):
        return self.variant_site_marker is not None

    def __str__(self):
        return '[{alleles}]'.format(alleles='|'.join(self.alleles))


def sites_mask(regions):
    """Generate sites masks."""
    sites = []

    for region in regions:
        if len(region.alleles) == 1:
            sites += ['0'] * len(region.alleles[0])
            continue

        sites.append('0')
        for i, allele in enumerate(region.alleles):
            sites += [str(region.variant_site_marker)] * len(allele)
            sites.append('0')

    return '\t'.join(sites)


def alleles_mask(regions):
    """Generate alleles masks."""
    alleles = []

    for region in regions:
        if len(region.alleles) == 1:
            alleles += ['0'] * len(region.alleles[0])
            continue

        alleles.append('0')
        for i, allele in enumerate(region.alleles):
            alleles += [str(i + 1)] * len(allele)
            alleles.append('0')

    return '\t'.join(alleles)
