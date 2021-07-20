#include <stdio.h>
#include "ising.h"
#include "record.h"
#include "parse_args.h"

int main(int argc, char **argv)
{
    if (argc != 6) {
        fprintf(stderr, "usage: %s <j> <beta> <iterations> <count> <filename>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    double j = parseDouble(argv[1], "j");
    double beta = parseDouble(argv[2], "beta");
    unsigned long iterations = parseUnsignedLong(argv[3], "iterations");
    unsigned long count = parseUnsignedLong(argv[4], "count");

    char *filename = argv[5];

    FILE *data_file = fopen(filename, "w");
    if (!data_file) {
        perror("error opening file");
        exit(EXIT_FAILURE);
    }
    
    state_t lattice[TIME_LEN * SPINS_PER_STATE_T];

    writeHeader(data_file, j, beta);
    for (unsigned long i = 0; i < count; i++) {
        initLattice(lattice);
        metropolis(lattice, hamiltonian(lattice, j), j, beta, iterations);
        writeState(data_file, lattice);
    }
    fclose(data_file);
}
