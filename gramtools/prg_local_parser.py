## @file
# Extracts characters one by one from the prg, keeping only chunks of characters in memory.
# The characters are directly comparison to an allele_index, which is used to choose which allele needs to be picked for each variant site.

PRG_READ_CHUNK_SIZE = 1000000
DIGITS = {
    '0', '1', '2', '3',
    '4', '5', '6', '7',
    '8', '9'
}
FASTA_LINE_SIZE = 60  #In characters


## An OOP implementation of cursor-based 'local' prg parsing.
# @param allele_indexes a **generator** giving the allele index to choose for each variant site of the prg in `prg_fpath`
class Prg_Local_Parser(object):
    def __init__(self, prg_fpath, output_file_name, fasta_header, allele_indexes):
        self.prg_parser = _parse(prg_fpath)
        self.writer = FastaWriter(output_file_name, description = fasta_header)
        self.allele_indexes = allele_indexes

    def parse(self):
        _dump_fasta(self.prg_parser, self.allele_indexes, self.writer)
        self.writer.close()


class FastaWriter:
    def __init__(self, fpath, description):
        self._fpath = fpath
        self._fhandle = open(self._fpath, 'w')
        self._fhandle.write("> {} \n".format(description))
        self._running_tally = 0

        self._cache = []
        self._max_cache_size = 1000000

    def append(self, char):
        self._cache.append(char)
        self._running_tally += 1

        if self._running_tally == FASTA_LINE_SIZE:
            self._cache.append('\n')
            self._running_tally = 0

        if len(self._cache) > self._max_cache_size:
            self._flush()

    def _flush(self):
        cache = ''.join(self._cache)
        self._fhandle.write(cache)
        self._cache = []

    def _flush_endFile(self):
        if self._cache[-1] != '\n':
            self._cache.append('\n')

        cache = ''.join(self._cache)
        self._fhandle.write(cache)

    def close(self):
        self._flush_endFile()
        self._fhandle.close()



def _is_int(data):
    if data is None or data == '':
        return False
    if len(data) > 1:
        return True
    return data in DIGITS


## A cursor on a single character in the prg.
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
            if not char:
                break
            extra_chars.append(char)
        chars += ''.join(extra_chars)
    return chars


def _parse(prg_fpath):
    cursor = _Cursor()
    with open(prg_fpath) as file_handle:
        while True:
            chars = _read_chunk(file_handle)
            if chars is None:
                break

            for cursor in _parse_prg_structure(chars, cursor):
                yield cursor

def _dump_fasta(prg_parser, allele_indexes, writer):
    allele_index = next(allele_indexes)

    for cursor in prg_parser:
        if cursor.just_left_site:
            try:
                allele_index = next(allele_indexes)
            except StopIteration:
                allele_index = 0

        if cursor.on_marker:
            continue

        if not cursor.within_allele:
            writer.append(cursor.char)
            continue

        if cursor.allele_id == allele_index:
            writer.append(cursor.char)

