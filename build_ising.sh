#!/bin/bash
cd npy_array
make
cd ..
gcc ising.c -L./npy_array -l:libnpy_array.a -lm -Wall
