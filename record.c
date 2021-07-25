#include "endian.h"
#include "record.h"
#include "ising.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

const int8_t *file_ver_identifier = "ISI\x01";

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

npy_array_t createNpyDoubleArrayNd(int ndim, ...)
{
    va_list vararg;
    va_start(vararg, ndim);

    int endianness = !((uint8_t *)(&(int){1}))[0];
    npy_array_t result = (npy_array_t) { .data = NULL,
        .ndim = ndim, .endianness = '<' + 2 * endianness, .typechar = 'f', 
        .elem_size = sizeof(double), .fortran_order = 0,
    };

    if (ndim >= sizeof(result.shape) / sizeof(size_t)) {
        fprintf(stderr, "too many dimensions in ndarray\n");
        exit(EXIT_FAILURE);
    }

    size_t total_count = 1;
    for (int i = 0; i < ndim; i++) {
        int dim = va_arg(vararg, int);
        if (dim <= 0) {
            fprintf(stderr, "error, negative or zero dimension\n");
            exit(EXIT_FAILURE);
        }
        if ((total_count * dim) / dim != total_count) {
            fprintf(stderr, "error, multiplication overflow when calculating size\n");
            exit(EXIT_FAILURE);
        }
        result.shape[i] = dim;
        total_count *= dim;
    }
    va_end(vararg);

    result.data = calloc(total_count, sizeof(double));
    if (!result.data) {
        fprintf(stderr, "Error allocating numpy array, abort!\n");
        exit(EXIT_FAILURE);
    }
    return result;
}

int writeHeader(FILE *fp, double j, double beta)
{
    if (fwrite(file_ver_identifier, 1, sizeof(file_ver_identifier) - 1, fp) != sizeof(file_ver_identifier) - 1)
        return -1;

    // make little endian versions of each to write to file
    // yes, this is ugly, but the alternatives are not as portable
    uint64_t j_le = htole64(*((uint64_t *)&j));
    if (fwrite(&j_le, sizeof(j_le), 1, fp) != 1)
        return -1;
    uint64_t beta_le = htole64(*((uint64_t *)&beta));
    if (fwrite(&beta_le, sizeof(beta_le), 1, fp) != 1)
        return -1;
    uint16_t lattice_time_len = htole16(TIME_LEN);
    if (fwrite(&lattice_time_len, sizeof(lattice_time_len), 1, fp) != 1)
        return -1;
    uint16_t lattice_space_len = htole16(SPACE_LEN);
    if (fwrite(&lattice_space_len, sizeof(lattice_space_len), 1, fp) != 1)
        return -1;
    uint16_t state_t_bytes = htole16(sizeof(state_t) * CHAR_BIT / 8);
    if (fwrite(&state_t_bytes, sizeof(state_t_bytes), 1, fp) != 1)
        return -1;

    return 0;
}

// writes a lattice to the specified FILE *
// using little-endian byte ordering
int writeState(FILE *fp, state_t *lattice)
{
    state_t buffer[SPACE_STATE_COUNT * TIME_LEN];
    for (size_t i = 0; i < SPACE_STATE_COUNT * TIME_LEN; i++) {
        // this will get optimized away by the compiler
        switch (sizeof(state_t) * CHAR_BIT) {
        case 64:
            buffer[i] = htole64(lattice[i]);
            break;
        case 32:
            buffer[i] = htole32(lattice[i]);
            break;
        case 16:
            buffer[i] = htole16(lattice[i]);
            break;
        default:
            buffer[i] = lattice[i];
            break;
        }
    }

    if (fwrite(buffer, sizeof(state_t), SPACE_STATE_COUNT * TIME_LEN, fp) != SPACE_STATE_COUNT * TIME_LEN)
        return -1;

    return 0;
}

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
    if (fread(&j_le, sizeof(j_le), 1, fp) != 1)
        return ERROR_READ;
    *(uint64_t *)j = le64toh(j_le);

    uint64_t beta_le;
    if (fread(&beta_le, sizeof(beta_le), 1, fp) != 1)
        return ERROR_READ;
    *(uint64_t *)beta = le64toh(beta_le);

    uint16_t lattice_time_len;
    if (fread(&lattice_time_len, sizeof(lattice_time_len), 1, fp) != 1)
        return ERROR_READ;
    if (lattice_time_len != htole16(TIME_LEN))
        return ERROR_LATTICE_SIZE;

    uint16_t lattice_space_len;
    if (fread(&lattice_space_len, sizeof(lattice_space_len), 1, fp) != 1)
        return ERROR_READ;
    if (lattice_space_len != htole16(SPACE_LEN))
        return ERROR_LATTICE_SIZE;

    // could probably adapt for mismatched size at runtime but I don't really feel like it
    uint16_t state_t_bytes; 
    if (fread(&state_t_bytes, sizeof(state_t_bytes), 1, fp) != 1)
        return ERROR_READ;
    if (state_t_bytes != htole16(sizeof(state_t) * CHAR_BIT / 8))
        return ERROR_STATE_SIZE;

    return READ_SUCCESS;
}

int readState(FILE *fp, state_t *lattice)
{
    if (fread(lattice, sizeof(state_t), SPACE_STATE_COUNT * TIME_LEN, fp) != SPACE_STATE_COUNT * TIME_LEN)
        return ERROR_READ;
    
    for (int i = 0; i < SPACE_STATE_COUNT * TIME_LEN; i++) {
        switch (sizeof(state_t) * CHAR_BIT) {
        case 64:
            lattice[i] = le64toh(lattice[i]);
            break;
        case 32:
            lattice[i] = le32toh(lattice[i]);
            break;
        case 16:
            lattice[i] = le16toh(lattice[i]);
            break;
        default:
            lattice[i] = lattice[i];
            break;
        }
    }

    return READ_SUCCESS;
}
