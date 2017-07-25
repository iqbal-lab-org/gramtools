class GenomeRegions:
    """Variantblocks container class. Contains entire genome
    split into variant blocks.
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
            idx_range = range(start_region_idx, len(self._regions))
        else:
            idx_range = range(start_region_idx, -1, -1)

        for idx in idx_range:
            region = self._regions[idx]
            yield region

    def add_region(self, region, current_var_marker):
        """Add a region of the genome."""
        region = _GenomeRegion(region, current_var_marker)
        self._regions.append(region)
        self._region_idx[region] = len(self._regions) - 1


class _GenomeRegion:
    """A region of the genome, either a variant site or a non-variant site."""

    def __init__(self, alleles, current_var_marker=None):
        # TODO: why is alleles copied here?
        self.alleles = list(alleles)
        self.min_allele_len = min(len(i) for i in alleles)
        self.variant_site_marker = current_var_marker

    @property
    def is_variant_site(self):
        return self.variant_site_marker is not None


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
