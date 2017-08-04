def _handle_var_region(region, prg, marker):
    prg += str(marker)
    for i, allele in enumerate(region):
        prg += allele
        at_last_allele = (i == len(region) - 1)
        if at_last_allele:
            break
        prg += str(marker + 1)
    prg += str(marker)
    marker += 2
    return prg, marker


def compose_prg(prg_structure):
    """Given a prg structure return a prg string.

    prg_structure = [
        ['AC'],
        ['TA', 'G'],
        ['C'],
        ['A', 'T'],
        ['GG'],
    ]

    Return value: AC5TA6G5C7A8T7GG
    """
    prg = ''
    marker = 5
    for region in prg_structure:
        is_var_region = (len(region) > 1)
        if is_var_region:
            prg, marker = _handle_var_region(region, prg, marker)
        else:
            prg += region[0]
    return prg