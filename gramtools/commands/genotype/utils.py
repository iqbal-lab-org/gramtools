import json


def _load_grouped_allele_coverage(fpath):
    with open(fpath, "r") as fhandle:
        data = json.load(fhandle)
    groups_coverage = data["grouped_allele_counts"]

    allele_groups = groups_coverage["allele_groups"]
    for key, value in allele_groups.items():
        allele_groups[key] = set(value)

    groups_site_counts = groups_coverage["site_counts"]
    return allele_groups, groups_site_counts


def _load_per_base_coverage(fpath):
    with open(fpath, "r") as fhandle:
        data = json.load(fhandle)
    data = data["allele_base_counts"]
    return data
