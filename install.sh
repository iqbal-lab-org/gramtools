#!/bin/bash

cd sdsl-lite
./install.sh ..

cd ..
g++ -isystem ./googletest/include -I./googletest -pthread -c ./googletest/src/gtest-all.cc
ar -rv ./lib/libgtest.a gtest-all.o