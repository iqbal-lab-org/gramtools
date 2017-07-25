from . import genome_regions


var_marker_chars = set('0123456789')


class CurrentState:
    """Current cursor state."""
    def __init__(self):
        self.in_var_site = False
        self.current_var_marker = None
        self.region = []


def flush_region(state, regions):
    """Reset state for next region."""
    if not state.region:
        return
    state.region = [''.join(x) for x in state.region]
    regions.add_region(state.region, state.current_var_marker)
    state.region = []


def handle_var_site(prg_char, state, regions):
    """Handle cursor in variant site."""
    is_marker_char = prg_char[0] in var_marker_chars

    entering_var_site = (not state.in_var_site
                         and is_marker_char)
    if entering_var_site:
        flush_region(state, regions)
        state.in_var_site = True
        state.current_var_marker = prg_char
        return

    leaving_var_site = prg_char == state.current_var_marker
    if leaving_var_site:
        flush_region(state, regions)
        state.in_var_site = False
        state.current_var_marker = None
        return

    entering_next_allele = is_marker_char and state.in_var_site
    if entering_next_allele:
        state.region.append([])
        return

    in_allele = not is_marker_char and state.in_var_site
    if in_allele:
        if not state.region:
            state.region.append([])
        state.region[-1].append(prg_char)
        return


def handle_nonvar_site(prg_char, state):
    """Handle cursor in non-variant site."""
    is_nonvar_base = (not state.in_var_site
                      and prg_char[0] not in var_marker_chars)
    if is_nonvar_base:
        if not state.region:
            state.region.append([])
        state.region[-1].append(prg_char)
        return


class IterPeek:
    """An iterator which stores the next (future) value for peeking."""

    def __init__(self, iterable):
        self.iter = iter(iterable)
        self.peek = None
        self.stop_flag = False

    def __iter__(self):
        return self

    def __next__(self):
        if self.stop_flag:
            raise StopIteration

        if self.peek is None:
            current = next(self.iter)
        else:
            current = self.peek

        try:
            self.peek = next(self.iter)
        except StopIteration:
            self.peek = None
            self.stop_flag = True
        return current


def decode_prg(prg):
    """Decode prg by concatenating variant site marker digits."""
    iter_prg = IterPeek(prg)
    marker = ''
    for x in iter_prg:
        if x not in var_marker_chars:
            yield x
            continue

        marker += x
        if iter_prg.peek in var_marker_chars:
            continue
        else:
            x = marker
            marker = ''
        yield x


def parse(prg):
    """Process genome _regions."""
    regions = genome_regions.GenomeRegions()
    state = CurrentState()
    for x in decode_prg(prg):
        handle_var_site(x, state, regions)
        handle_nonvar_site(x, state)

    flush_region(state, regions)
    return regions
