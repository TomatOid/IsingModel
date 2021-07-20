#pragma once
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <endian.h>
#include "npy_array/npy_array.h"

npy_array_t createNpyDoubleArray1D(size_t count)
{
    // on little endian, this returns 0, on big endian, this is 1
    // this will compile error out on super weird machines that have 
    // sizeof(char) == sizeof(int) instead of returning false results
    int endianness = !((uint8_t *)(&(int){1}))[0];
    npy_array_t result = (npy_array_t) { .data = calloc(count, sizeof(double)),
        .ndim = 1, .endianness = '<' + 2 * endianness, .typechar = 'f', 
        .elem_size = sizeof(double), .fortran_order = 0,
    };
    if (!result.data) {
        fprintf(stderr, "Error allocating numpy array, abort!\n");
        exit(EXIT_FAILURE);
    }
    result.shape[0] = count;
    return result;
}

const uint8_t *file_ver_identifier = "ISI\x01";

#define FILE_HEADER_LEN 26

int writeHeader(FILE *fp, double j, double beta)
{
    if (fwrite(file_ver_identifier, 1, sizeof(file_ver_identifier) - 1, fp) != sizeof(file_ver_identifier) - 1)
        return -1;

    // make little endian versions of each to write to file
    // yes, this is ugly, but the alternatives are not as portable
    uint64_t j_le = htole64(j);
    if (fwrite(&j_le, sizeof(j_le), 1, fp) != sizeof(j_le))
        return -1;
    uint64_t beta_le = htole64(beta);
    if (fwrite(&beta_le, sizeof(beta_le), 1, fp) != sizeof(beta_le))
        return -1;
    uint16_t lattice_time_len = htole16(TIME_LEN);
    if (fwrite(&lattice_time_len, sizeof(lattice_time_len), 1, fp) != sizeof(lattice_time_len))
        return -1;
    uint16_t lattice_space_len = htole16(SPACE_LEN);
    if (fwrite(&lattice_space_len, sizeof(lattice_space_len), 1, fp) != sizeof(lattice_space_len))
        return -1;
    uint16_t state_t_bytes = htole16(sizeof(state_t) * CHAR_BIT / 8);
    if (fwrite(&state_t_bytes, sizeof(state_t_bytes), 1, fp) != sizeof(state_t_bytes))
        return -1;

    return 0;
}

// writes a lattice to the specified FILE *
// using little-endian byte ordering
int writeState(FILE *fp, state_t *lattice)
{
    char buffer[SPACE_STATE_COUNT * TIME_LEN * sizeof(state_t)];
    for (size_t i = 0; i < SPACE_STATE_COUNT * TIME_LEN; i++) {
#if (sizeof(state_t) * CHAR_BIT) == 64
        buffer[i] = htole64(lattice[i]);
#elif (sizeof(state_t) * CHAR_BIT) == 32
        buffer[i] = htole32(lattice[i]);
#elif (sizeof(state_t) * CHAR_BIT) == 16
        buffer[i] = htole16(lattice[i]);
#else
        buffer[i] = lattice[i];
#endif
    }

    if (fwrite(buffer, sizeof(state_t), SPACE_STATE_COUNT * TIME_LEN, fp) != SPACE_STATE_COUNT * TIME_LEN)
        return -1;

    return 0;
}

enum {
    READ_SUCCESS,
    ERROR_BAD_PREFIX,
    ERROR_LATTICE_SIZE,
    ERROR_STATE_SIZE,
    ERROR_READ,
};

// read the lattice from the specified FILE pointer
// assumes all pointers are valid, returns 0 on success
int readHeader(FILE *fp, double *j, double *beta)
{
    uint8_t id_str[sizeof(file_ver_identifier) - 1];
    if (fread(id_str, 1, sizeof(id_str), fp) != sizeof(id_str))
        return ERROR_READ;

    if (memcmp(id_str, file_ver_identifier, sizeof(id_str)) != 0)
        return ERROR_BAD_PREFIX;

    uint64_t j_le;
    if (fread(&j_le, sizeof(j_le), 1, fp) != sizeof(j_le))
        return ERROR_READ;
    *j = (double)le64toh(j_le);

    uint64_t beta_le;
    if (fread(&beta_le, sizeof(beta_le), 1, fp) != sizeof(beta_le))
        return ERROR_READ;
    *beta = (double)le64toh(beta_le);

    uint16_t lattice_time_len;
    if (fread(&lattice_time_len, sizeof(lattice_time_len), 1, fp) != sizeof(lattice_time_len))
        return ERROR_READ;
    if (lattice_time_len != htole16(TIME_LEN))
        return ERROR_LATTICE_SIZE;

    uint16_t lattice_space_len;
    if (fread(&lattice_space_len, sizeof(lattice_space_len), 1, fp) != sizeof(lattice_space_len))
        return ERROR_READ;
    if (lattice_space_len != htole16(SPACE_LEN))
        return ERROR_LATTICE_SIZE;

    // could probably adapt for mismatched size at runtime but I don't really feel like it
    uint16_t state_t_bytes; 
    if (fread(&state_t_bytes, sizeof(state_t_bytes), 1, fp) != sizeof(state_t_bytes))
        return ERROR_READ;
    if (state_t_bytes != htole16(sizeof(state_t) * CHAR_BIT / 8))
        return ERROR_STATE_SIZE;

    return READ_SUCCESS;
}

int readLattice(FILE *fp, state_t *lattice)
{
    if (fread(lattice, sizeof(state_t), SPACE_STATE_COUNT * TIME_LEN, fp) != SPACE_STATE_COUNT * TIME_LEN)
        return ERROR_READ;
    
    for (int i = 0; i < SPACE_STATE_COUNT * TIME_LEN; i++) {
#if (sizeof(state_t) * CHAR_BIT) == 64
        lattice[i] = le64toh(lattice[i]);
#elif (sizeof(state_t) * CHAR_BIT) == 32
        lattice[i] = le32toh(lattice[i]);
#elif (sizeof(state_t) * CHAR_BIT) == 16
        lattice[i] = le16toh(lattice[i]);
#else
        lattice[i] = lattice[i];
#endif
    }

    return READ_SUCCESS;
}
