#include "ising.h"
#include "npy_array/npy_array.h"
#include "record.h"
#include "parse_args.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

int main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr, "usage: %s <j> <h*mu> <beta> <iterations>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    double j    = parseDouble(argv[1], "j");
    double h_mu = parseDouble(argv[2], "h_mu");
    double beta = parseDouble(argv[3], "beta");
    
    unsigned long iterations = parseUnsignedLong(argv[4], "iterations");

    state_t hot_lattice[TIME_LEN * SPACE_STATE_COUNT]  = { 0 };
    state_t cold_lattice[TIME_LEN * SPACE_STATE_COUNT] = { 0 };

    // fill hot_lattice with random spins
    initLattice(hot_lattice);

    npy_array_t hot_energies  = createNpyDoubleArray1D(iterations);
    npy_array_t cold_energies = createNpyDoubleArray1D(iterations);

    // iterate the metropolis algorithm, saving the energies of the hot and cold
    // lattices to the numpy arrays for later graphing
    double hot_energy  = hamiltonian(hot_lattice, j, h_mu);
    double cold_energy = hamiltonian(cold_lattice, j, h_mu);
    for (unsigned long i = 0; i < iterations; i++) {
        ((double *)(hot_energies.data))[i] = hot_energy / (SPACE_LEN * TIME_LEN);
        hot_energy  = metropolis(hot_lattice, hot_energy, j, h_mu, beta, SPACE_LEN * TIME_LEN);
        ((double *)(cold_energies.data))[i] = cold_energy / (SPACE_LEN * TIME_LEN);
        cold_energy = metropolis(cold_lattice, cold_energy, j, h_mu, beta, SPACE_LEN * TIME_LEN);
    }

    npy_array_save("hot_energies.npy", &hot_energies);
    npy_array_save("cold_energies.npy", &cold_energies);

    printf("hot: %f, cold: %f\n", hot_energy, cold_energy);

    return 0;
}
