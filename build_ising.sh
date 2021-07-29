#!/bin/bash
cd npy_array
make
cd ..
gcc -c ising.c record.c -fopenmp -g -O3
gcc ising.o record.o hot_v_cold.c -L./npy_array -l:libnpy_array.a -lm -Wall -o hot_v_cold
gcc ising.o record.o generate_states.c -O3 -fopenmp -lm -Wall -o generate_states
gcc ising.o record.o correlation.c -g -O3 -L./npy_array -l:libnpy_array.a -lm -o correlation
