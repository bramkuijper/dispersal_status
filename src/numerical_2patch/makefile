# makefile to compile the two-patch simulation

EXE=twopatch.exe
CPP=twopatch.cpp solve.cpp
HPP=twopatch.hpp
CFLAGS=-Wall -ggdb -std=c++17
CXX=g++
LIBS=-lgsl -lgslcblas


$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CFLAGS) -o twopatch.exe $(CPP) $(LIBS)

