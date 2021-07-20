#include "ising.h"
#include "npy_array/npy_array.h"
#include "record.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr, "usage: %s <j> <beta> <iterations>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *end_ptr;
    double j = strtod(argv[1], &end_ptr);
    if (errno == ERANGE || end_ptr == argv[1]) {
        fprintf(stderr, "j must be a valid double-width floating point number, got %s\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    double beta = strtod(argv[2], &end_ptr);
    if (errno == ERANGE || end_ptr == argv[2]) {
        fprintf(stderr, "beta must be a valid double-width floating point number, got %s\n", argv[2]);
        exit(EXIT_FAILURE);
    }
    unsigned long iterations = strtoul(argv[3], &end_ptr, 10);
    if (errno == ERANGE) {
        fprintf(stderr, "value %s for iterations is out of range for type unsigned long\n", argv[3]);
        exit(EXIT_FAILURE);
    }
    else if (errno == EINVAL || end_ptr == argv[3]) {
        fprintf(stderr, "iterations must be a valid positive integer, got %s\n", argv[3]);
        exit(EXIT_FAILURE);
    }

    state_t hot_lattice[TIME_LEN * SPACE_STATE_COUNT]  = { 0 };
    state_t cold_lattice[TIME_LEN * SPACE_STATE_COUNT] = { 0 };

    // fill hot_lattice with random spins
    initLattice(hot_lattice);

    npy_array_t hot_energies  = createNpyDoubleArray1D(iterations);
    npy_array_t cold_energies = createNpyDoubleArray1D(iterations);

    // iterate the metropolis algorithm, saving the energies of the hot and cold
    // lattices to the numpy arrays for later graphing
    double hot_energy  = hamiltonian(hot_lattice, j);
    double cold_energy = hamiltonian(cold_lattice, j);
    for (unsigned long i = 0; i < iterations; i++) {
        ((double *)(hot_energies.data))[i] = hot_energy / (SPACE_LEN * TIME_LEN);
        hot_energy  = metropolis(hot_lattice, hot_energy, j, beta, SPACE_LEN * TIME_LEN);
        ((double *)(cold_energies.data))[i] = cold_energy / (SPACE_LEN * TIME_LEN);
        cold_energy = metropolis(cold_lattice, cold_energy, j, beta, SPACE_LEN * TIME_LEN);
    }

    npy_array_save("hot_energies.npy", &hot_energies);
    npy_array_save("cold_energies.npy", &cold_energies);

    printf("hot: %f, cold: %f\n", hot_energy, cold_energy);

    return 0;
}
