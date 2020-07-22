CXX=g++
CXXFLAGS=-Wall -O3 -std=c++17
EXE1=disp_status_2p.exe
EXE_TEST=testcases.exe

all : $(EXE1) $(EXE_TEST)

$(EXE1) : dispersal_status_2p.hpp dispersal_status_2p.cpp
	$(CXX) $(CXXFLAGS) -o $(EXE) dispersal_status_2p.cpp

$(EXE_TEST) : unittesting.cpp patchfreq.cpp patchfreq.hpp catch.hpp
	$(CXX) $(CXXFLAGS) -o $(EXE_TEST) unittesting.cpp patchfreq.cpp

clean :
	rm -rf $(EXE)
	rm -rf $(EXE_TEST)
