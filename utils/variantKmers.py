import re
from Bio import SeqIO
import argparse
import array
import sys

class Chain:
        '''Chain of Variantblocks, container  class, contains all genome split into VBlocks, like a linked list'''
	def __init__(self):
		self.blcks={}          # Map to address each block by id (id are odd numbers)
		self.firstBlock=None   # Pointer to first block
		self.lastBlock=None    # Pointer to last block

	def addVariantBlock(self,b):   # Adding  a block to the linked list
	    if not self.firstBlock:
		self.firstBlock=self.lastBlock=b
	    else:
		self.lastBlock.next=b
		b.prev=self.lastBlock
		self.lastBlock=b

	    if b.bid!=0:
		self.blcks[b.bid]=b

	def __getitem__(self,bid):       # Getting a Vblock by id
	    return self.blcks[bid]

	def __iter__(self):             # Iterating over all blocks
	    a=self.firstBlock
	    while a:
		yield a
		a=a.next

	def blocksAtKDistance(self,bl,k):  
            '''given a distance (k:kmer-length) which is the further block in each direction from a given block?
               this method answers that question'''
	    first=bl
	    dist=k-1
	    while dist>0 and first.prev:
		first=first.prev
		dist-=len(first)

	    last=bl
	    dist=k-1
	    while dist>0 and last.next:
		last=last.next
		dist-=len(last)

	    return first,last

	def stringsForVariant(self,bid,k): 
            '''Given a block, and the two farther blocks in each direction, computes all possible strings in volving variants
                it's a generator, the strings are not generated all at once! :D '''
	    block=self[bid]
	    fb,lb=self.blocksAtKDistance(block,k)
	    if lb: lb=lb.next

	    stringsStart=set()
	    stringsEnd=set()

	    for i in fb.getStringsUntil(block):
		stringsStart.add(i[-k+1:].tostring())
	    if block.next: 
		for i in block.next.getStringsUntil(lb):
		    stringsEnd.add(i[:k-1].tostring())

	    for s in stringsStart:
		for c in block:
		    for e in stringsEnd:
			yield s+c.tostring()+e

	def kmersForVariant(self,bid,k):
                '''Given a string, retirns all it's kmers'''
		for s in self.stringsForVariant(bid,k):
			for i in xrange (0,len(s)-k+1):
				yield s[i:i+k]


	def getKmerDict(self,k,nonVariants=False):
            '''Returns a dictionary containing a map kmer->set(variants) for all variants''' 
            s=set()
	    for i in self.blcks:
		for j in self.kmersForVariant(i,k):
                    if j not in s:
                       yield j
                       s.add(j)

            if nonVariants:
                q=self.firstBlock
                while q:
                    if not q.bid:
                        c=q.cads[0]
                        for i in xrange (0,len(c)-k+1):
                            j=c[i:i+k].tostring()
                            if j not in s:
                                yield j
                                s.add(j)
                    q=q.next

	 #   return res

        def parseString(self,cad):
                '''Feeds that chain with a fasta String'''
                start=0
                out=True
            #    for i in re.finditer("[^0-9][0-9]*[13579][^0-9]",cad):
            #    	if out:
            #    		self.addVariantBlock( VariantBlock([cad[start:i.start()+1]]))
            #    		out=False
            #    		start=i.start()+1
            #    	else:
            #    		endpos=i.span()[1]-1
            #    		self.addVariantBlock(VariantBlock(re.findall("[^0-9]+",cad[start:endpos]),int(cad[i.span()[0]+1:i.span()[1]-1])))
            #    		start=endpos
            #    		out=True
            #    self.addVariantBlock(VariantBlock([cad[start:]]))
                i=0
                lastblock=0
                while i<len(cad):
                    if cad[i] not in "0123456789": 
                        i+=1
                        continue
                    lastsequence=cad[lastblock:i]
                    if lastsequence: self.addVariantBlock(VariantBlock([lastsequence],0))

                    ns=""
                    blockini=i
                    while cad[i] in "0123456789":
                        ns+=cad[i]
                        i+=1
                    # We look for the end of the block (end of the marker)
                    nextn=cad.find(ns,i)
                    # we jump the marker (the odd number)
                    lastblock=i=nextn+len(ns)
                    blockend=i
                    self.addVariantBlock(VariantBlock(re.findall("[^0-9]+",cad[blockini:blockend]),int(ns)))
                        

                lastsequence=cad[lastblock:]
                if lastsequence: self.addVariantBlock(VariantBlock([lastsequence],0))


            

                

        def parseFasta(self,f):
                '''Feeds that chain with a fasta file'''
                return self.parseString(str(SeqIO.read(f,'fasta').seq))


class VariantBlock:
	def __init__(self,cads=[],bid=0):
		self.cads=[array.array('c',i) for i in cads]
		self.next=None
		self.prev=None
		self.bid=bid
		self.lgth=min([len(i) for i in cads])


	def __len__(self):
		return self.lgth

	def getKmer(self,pos,lgt):
		for i in cads:
			if len(i): pass

	def __iter__(self):
		for i in self.cads:
		    yield i

	def __str__(self):
		return "{0} ({2}): {1}".format(self.bid,self.cads,len(self))

	def getStringsUntil(self,end):
		for i in self.cads:
			if self.next and self.next!=end:
				for j in self.next.getStringsUntil(end):
					yield i+j
			else:
				yield i

	__repr__=__str__



if __name__=="__main__":
        aparser = argparse.ArgumentParser(description='PRG kmer generator') 
        aparser.add_argument("-f", dest="fasta",  help="Fasta file", required=True) 
        aparser.add_argument("-n", dest="nonvariant", default=False, help="Print kmers from intervariant regions", action='store_true') 
#        aparser.add_argument("-v", dest="vnumber",help="output variant numbers associated",default=False, action='store_true') 
        aparser.add_argument("-k", dest="ksize",help="kmer size [31]", type=int, default=31) 
        options = aparser.parse_args()

        bls=Chain()


        bls.parseFasta(options.fasta)
        
        for i in  bls.getKmerDict(options.ksize,options.nonvariant):
            print i
