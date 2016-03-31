CC=g++

CPPFLAGS=-std=c++11 -g 
INCL=-I./htslib -I. 
LDFLAGS=-lz


default: seq

gzstream.o: gzstream.C gzstream.h
	${CXX} -c -o gzstream.o gzstream.C $(INCL) 

seq: seqread.cpp gzstream.o
	$(CC) -o seq seqread.cpp gzstream.o  $(INCL) $(LDFLAGS) $(CPPFLAGS)


clean:
	rm seq *.o
