CC=g++

CPPFLAGS=-std=c++11 -g 
INCL=-I./htslib -I. 
LDFLAGS=-lz
HTSINC=-I./htslib -L./htslib -lz -lpthread


default: seq

gzstream.o: gzstream.C gzstream.h
	${CXX} -c -o gzstream.o gzstream.C $(INCL) 

seq: seqread.cpp gzstream.o
	$(CC) -o seq seqread.cpp gzstream.o  $(INCL) $(LDFLAGS) $(CPPFLAGS)

test: test.c
	$(CC) $(HTSINC) test.c -o test

clean:
	rm seq *.o
