import re
import argparse
import array
import random

from Bio import SeqIO


class Chain:
    """Chain of Variantblocks, container class, contains all genome split
    into VBlocks, like a linked list.
    """
    def __init__(self):
        # Map to address each block by id (id are odd numbers)
        self.blcks = {}
        # Pointer to first block
        self.first_block = None
        # Pointer to last block
        self.last_block = None

    def add_variant_block(self, block):
        """Add block to linked list."""
        if not self.first_block:
            self.first_block = self.last_block = block
        else:
            self.last_block.next = block
            block.prev = self.last_block
            self.last_block = block

        if block.bid != 0:
            self.blcks[block.bid] = block

    def __getitem__(self, bid):
        """Getting a Vblock by id."""
        return self.blcks[bid]

    def __iter__(self):
        """Iterating over all blocks."""
        block = self.first_block
        while block:
            yield block
            block = block.next

    @staticmethod
    def blocks_at_kdistance(bl, k):
        """Given a distance (k:kmer-length) which is the further
        block in each direction from a given block? this method answers
        that question.
        """
        first = bl
        dist = k-1
        while dist > 0 and first.prev:
            first = first.prev
            dist -= len(first)

        last = bl
        dist = k-1
        while dist > 0 and last.next:
            last = last.next
            dist -= len(last)

        return first, last

    def strings_for_variant(self, bid, k):
        """Given a block, and the two farther blocks in each direction,
        computes all possible strings in volving variants
        it's a generator, the strings are not generated all at once!
        """
        block = self[bid]
        fb, lb = self.blocks_at_kdistance(block, k)
        if lb:
            lb = lb.next

        strings_start = set()
        strings_end = set()

        for string in fb.get_strings_until(block):
            strings_start.add(string[-k+1:].tostring())
        if block.next:
            for string in block.next.get_strings_until(lb):
                strings_end.add(string[:k-1].tostring())

        for s in strings_start:
            for c in block:
                for e in strings_end:
                    yield s+c.tostring()+e

    def kmers_for_variant(self, bid, k):
        """Given a string, retirns all it's kmers."""
        for s in self.strings_for_variant(bid, k):
            for x in xrange(0, len(s)-k+1):
                yield s[x: x + k]

    def random_genome(self):
        genome = []
        f = self.first_block
        while f:
            genome.append(random.choice(f.cads).tostring())
            f = f.next
        return ''.join(genome)

    def get_mask(self):
        sites_mask = []
        alleles_mask = []
        f = self.first_block
        while f:
            if len(f.cads) == 1:
                sites_mask += ["0"]*len(f.cads[0])
                alleles_mask += ["0"]*len(f.cads[0])
            else:
                sites_mask.append("0")
                alleles_mask.append("0")
                k = 1
                for cad in f.cads:
                    sites_mask += [str(f.bid)]*len(cad)
                    alleles_mask += [str(k)]*len(cad)
                    sites_mask.append("0")
                    alleles_mask.append("0")
                    k += 1

            f = f.next
        return "\t".join(sites_mask), "\t".join(alleles_mask)

    def get_kmer_dict(self, k, non_variants=False):
        """Return a dictionary containing a map kmer->set(variants)
        for all variants.
        """
        s = set()
        for block in self.blcks:
            for j in self.kmers_for_variant(block, k):
                if j not in s:
                    yield j
                    s.add(j)

                if non_variants:
                    q = self.first_block
                    while q:
                        if not q.bid:
                            c = q.cads[0]
                            for x in xrange(0, len(c) - k + 1):
                                j = c[x: x + k].tostring()
                                if j not in s:
                                    yield j
                                    s.add(j)
                        q = q.next

    def parse_string(self, cad):
        """Feeds that chain with a fasta String."""
        x = 0
        lastblock = 0
        while x < len(cad):
            if not cad[x].isdigit():
                x += 1
                continue
            lastsequence = cad[lastblock: x]
            if lastsequence:
                self.add_variant_block(VariantBlock([lastsequence], 0))

            ns = ""
            blockini = x
            while cad[x].isdigit():
                ns += cad[x]
                x += 1
            # We look for the end of the block (end of the marker)
            nextn = cad.find(ns, x)
            # we jump the marker (the odd number)
            lastblock = x = nextn+len(ns)
            blockend = x
            self.add_variant_block(VariantBlock(
                re.findall("[^0-9]+", cad[blockini: blockend]), int(ns)))

        lastsequence = cad[lastblock:]
        if lastsequence:
            self.add_variant_block(VariantBlock([lastsequence], 0))

    def parse_fasta(self, f):
        """Feeds that chain with a fasta file."""
        return self.parse_string(str(SeqIO.read(f, 'fasta').seq))


class VariantBlock:
    def __init__(self, cads=None, bid=0):
        if cads is None:
            cads = []
        self.cads = [array.array('c', cad) for cad in cads]
        self.lgth = min([len(cad) for cad in cads])
        self.next = None
        self.prev = None
        self.bid = bid

    def __len__(self):
        return self.lgth

    def get_kmer(self):
        for cad in self.cads:
            if len(cad):
                pass

    def __iter__(self):
        for cad in self.cads:
            yield cad

    def __str__(self):
        return "{0} ({2}): {1}".format(self.bid, self.cads, len(self))

    def get_strings_until(self, end):
        for cad in self.cads:
            if self.next and self.next != end:
                for j in self.next.get_strings_until(end):
                    yield cad + j
            else:
                yield cad

    __repr__ = __str__


if __name__ == "__main__":
    aparser = argparse.ArgumentParser(description='PRG kmer generator')
    aparser.add_argument("-f", dest="fasta", help="Fasta file",
                         required=True)
    aparser.add_argument("-n", dest="nonvariant", default=False,
                         help="Print kmers from intervariant regions",
                         action='store_true')
    aparser.add_argument("-k", dest="ksize", help="kmer size [31]",
                         type=int, default=31)
    aparser.add_argument("-r", dest="nreads",
                         help="Generates random genome and random reads",
                         type=int, default=0)
    aparser.add_argument("-m", dest="mask", action='store_true',
                         default=False, help="Dump mask (mask.txt)")
    options = aparser.parse_args()

    bls = Chain()
    bls.parse_fasta(options.fasta)

    for i in bls.get_kmer_dict(options.ksize, options.nonvariant):
        print(i)

    if options.nreads:
        rgen = bls.random_genome()
        a = open("randomGenome.fa", "w")
        a.write(">randGenome\n{0}\n".format(rgen))
        a.close()

        a = open("randomReads.fastq", "w")
        for i in xrange(options.nreads):
            posic = random.randint(0, len(rgen)-150)
            if random.randint(0, 1):
                # a.write("@{0}\n{1}\n+\n{2}\n".format(i,str(Seq.Seq(rgen[posic:posic+150]).reverse_complement()),'H'*150))
                a.write("@{0}\n{1}\n+\n{2}\n".format(
                    i, rgen[posic:posic+150], 'H'*150))
            else:
                a.write("@{0}\n{1}\n+\n{2}\n".format(
                    i, rgen[posic:posic+150], 'H'*150))

        a.close()

    if options.mask:
        sites, alleles = bls.get_mask()
        a = open('mask_sites.txt', 'w')
        a.write(sites)
        a.close()
        a = open('mask_alleles.txt', 'w')
        a.write(alleles)
        a.close()
