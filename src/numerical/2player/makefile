CXX=g++
CXXFLAGS=-Wall -O3 -std=c++17
EXE1=disp_status_2p.exe
EXE_TEST=testcases.exe

all : $(EXE1) $(EXE_TEST)

$(EXE1) : dispersal_status_2p.hpp dispersal_status_2p.cpp patchfreq.cpp patchfreq.hpp sol.cpp
	$(CXX) $(CXXFLAGS) -o $(EXE1) dispersal_status_2p.cpp patchfreq.cpp sol.cpp

$(EXE_TEST) : unittesting.cpp patchfreq.cpp patchfreq.hpp catch.hpp
	$(CXX) $(CXXFLAGS) -o $(EXE_TEST) unittesting.cpp patchfreq.cpp

clean :
	rm -rf $(EXE1)
	rm -rf $(EXE_TEST)
