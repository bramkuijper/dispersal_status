#!/usr/bin/env bash

./generate_cpp.py > numsolve2.cpp
g++ -Wall -o numsolve numsolve2.cpp -lm -lrt -lgsl -lgslcblas
