## @file
# Extracts characters one by one from the prg, keeping only chunks of characters in memory.
# The characters are amenable to direct use.

PRG_READ_CHUNK_SIZE = 10000
DIGITS = {
    '0', '1', '2', '3',
    '4', '5', '6', '7',
    '8', '9'
}


class FastaWriter:
    def __init__(self, fpath):
        self._fpath = fpath
        self._fhandle = open(self._fpath, 'w')
        self._fhandle.write('> A personal reference sequence built using gramtools.\n')

        self._cache = []
        self._max_cache_size = 10000

    def append(self, char):
        self._cache.append(char)
        if len(self._cache) > self._max_cache_size:
            self._flush()

    def _flush(self):
        cache = ''.join(self._cache)
        self._fhandle.write(cache)
        self._cache = []

    def close(self):
        self._flush()
        self._fhandle.close()



def _is_int(data):
    if data is None or data == '':
        return False
    if len(data) > 1:
        return True
    return data in DIGITS


def _parse_prg_chars(chars):
    int_chars = []
    for char in chars:
        if _is_int(char):
            int_chars.append(char)
            continue

        else:
            if int_chars:
                yield ''.join(int_chars)
                int_chars = []
            yield char

    if int_chars:
        yield ''.join(int_chars)


class _Cursor:
    char = None
    site_marker = None
    allele_id_counter = None

    on_marker = False
    just_left_site = False

    @property
    def allele_id(self):
        if not self.within_allele:
            return None
        return self.allele_id_counter

    @property
    def within_allele(self):
        return self.site_marker is not None and self.on_marker is False


def _parse_prg_structure(chars, cursor):
    for char in _parse_prg_chars(chars):
        cursor.char = char

        if cursor.just_left_site:
            cursor.just_left_site = False
            cursor.site_marker = None
            cursor.allele_id_counter = None

        if not _is_int(char):
            cursor.on_marker = False
            yield cursor
            continue

        cursor.on_marker = True

        site_boundary = int(char) % 2 != 0
        if site_boundary:
            entring = cursor.site_marker is None
            if entring:
                cursor.site_marker = int(char)
                cursor.allele_id_counter = 0
                yield cursor
                continue

            exiting = cursor.site_marker is not None
            if exiting:
                cursor.just_left_site = True
                cursor.allele_id_counter = None
                yield cursor
                continue

        allele_boundary = int(char) % 2 == 0
        if allele_boundary:
            cursor.allele_id_counter += 1
            yield cursor


def _read_chunk(file_handle, chunk_size=PRG_READ_CHUNK_SIZE):
    chars = file_handle.read(chunk_size)
    if not chars:
        return None

    if _is_int(chars[-1]):
        extra_chars = []
        while not extra_chars or _is_int(extra_chars[-1]):
            char = file_handle.read(1)
            if not chars:
                break
            extra_chars.append(char)
        chars += ''.join(extra_chars)
    return chars


def parse(prg_fpath):
    cursor = _Cursor()
    with open(prg_fpath) as file_handle:
        while True:
            chars = _read_chunk(file_handle)
            if chars is None:
                break

            for cursor in _parse_prg_structure(chars, cursor):
                yield cursor
