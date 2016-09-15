
# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = ./googletest

# Where to find user code.
USER_DIR = ./src
UTILS_INC=-I ./utils

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include $(UTILS_INC)

# Flags passed to the C++ compiler.
CXXFLAGS += -std=c++11 -g -Wall -Wextra -pthread 
EXTRAFLAGS += -O3

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = slow_unittest_bidir_search_bwd_fwd unittest_bidir_search_bwd_fwd 

OBJECTS=bidir_search.o skip.o get_location.o bidir_search_bwd.o \
	bidir_search_fwd.o precalc_kmer_matches.o parse_masks.o \
	csa_construction.o map.o 

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
#SDSL_HEADERS= ./include/
VBWT_HEADERS= ./include/
# If Boost is not installed system wide, specify the prefix here.
BOOST_HEADERS= /data2/apps/boost_1_60_0/

LIBS=./lib/ 

# House-keeping build targets.

all : $(TESTS) gramtools

clean :
	rm -f $(TESTS) gramtools *.o

%.o : $(USER_DIR)/%.cpp $(VBWT_HEADERS)/bwt_search.h 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXTRAFLAGS) -I $(BOOST_HEADERS) -I $(VBWT_HEADERS) \
		-L $(LIBS) -c $< -lsdsl -ldivsufsort -ldivsufsort64 -o $@

gramtools: $(OBJECTS) ./src/precalc_gen.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXTRAFLAGS) -I $(BOOST_HEADERS) -I $(VBWT_HEADERS) \
		-I $(BOOST_HEADERS) -L $(LIBS) $^ -o $@ \
		-lsdsl -ldivsufsort -ldivsufsort64 -lhts -lz

# TODO refactor these rules to make them more comprehensible.

unittest_bidir_search_bwd_fwd.o : $(USER_DIR)/test/unittest_bidir_search_bwd_fwd.cpp $(VBWT_HEADERS)/bwt_search.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I $(VBWT_HEADERS) -I $(BOOST_HEADERS) -c $(USER_DIR)/test/unittest_bidir_search_bwd_fwd.cpp

slow_unittest_bidir_search_bwd_fwd.o : $(USER_DIR)/test/slow_unittest_bidir_search_bwd_fwd.cpp $(VBWT_HEADERS)/bwt_search.h $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I $(VBWT_HEADERS) -I $(BOOST_HEADERS) -c $(USER_DIR)/test/slow_unittest_bidir_search_bwd_fwd.cpp

unittest_bidir_search_bwd_fwd : bidir_search.o skip.o get_location.o csa_construction.o bidir_search_bwd.o bidir_search_fwd.o unittest_bidir_search_bwd_fwd.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXTRAFLAGS) -I $(VBWT_HEADERS) -L $(LIBS) $^ -o $@ -lsdsl -ldivsufsort -ldivsufsort64 -lgtest -lpthread

slow_unittest_bidir_search_bwd_fwd : bidir_search.o skip.o get_location.o csa_construction.o bidir_search_bwd.o bidir_search_fwd.o slow_unittest_bidir_search_bwd_fwd.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXTRAFLAGS) -I $(VBWT_HEADERS) -L $(LIBS) $^ -o $@ -lsdsl -ldivsufsort -ldivsufsort64 -lgtest -lpthread
