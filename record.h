#pragma once
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include "ising.h"
#include "npy_array/npy_array.h"

#define FILE_HEADER_LEN 26

enum {
    READ_SUCCESS,
    ERROR_BAD_PREFIX,
    ERROR_LATTICE_SIZE,
    ERROR_STATE_SIZE,
    ERROR_READ,
};

npy_array_t createNpyDoubleArray1D(size_t count);

npy_array_t createNpyDoubleArrayNd(int count, ...);

npy_array_t createNpyArrayNd(char typechar, int type_size, int ndim, ...);

int writeHeader(FILE *fp, double j, double beta);

// writes a lattice to the specified FILE *
// using little-endian byte ordering
int writeState(FILE *fp, state_t *lattice);

// read the lattice from the specified FILE pointer
// assumes all pointers are valid, returns 0 on success
int readHeader(FILE *fp, double *j, double *beta);

int readState(FILE *fp, state_t *lattice);
