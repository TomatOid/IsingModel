#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "ising.h"
#include "record.h"
#include "parse_args.h"

// yoinked from David Blackman and Sebastiano Vigna's excellent 'PRNG shootout' page (CC0)
// equivalent to calling xorshift256() 2^128 times
void jump()
{
    static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    for(int i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
        for(int b = 0; b < 64; b++) {
            if (JUMP[i] & UINT64_C(1) << b) {
                s0 ^= xorshift_state[0];
                s1 ^= xorshift_state[1];
                s2 ^= xorshift_state[2];
                s3 ^= xorshift_state[3];
            }
            xorshift256();    
        }
        
    xorshift_state[0] = s0;
    xorshift_state[1] = s1;
    xorshift_state[2] = s2;
    xorshift_state[3] = s3;
}

int main(int argc, char **argv)
{
    if (argc != 7) {
        fprintf(stderr, "usage: %s <j> <h*mu> <beta> <iterations> <count> <filename>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    double j    = parseDouble(argv[1], "j");
    double h_mu = parseDouble(argv[2], "h_mu");
    double beta = parseDouble(argv[3], "beta");

    unsigned long iterations = parseUnsignedLong(argv[4], "iterations");
    unsigned long count      = parseUnsignedLong(argv[5], "count");

    char *filename = argv[6];

    FILE *data_file = fopen(filename, "w");
    if (!data_file) {
        perror("error opening file");
        exit(EXIT_FAILURE);
    }

    if (writeHeader(data_file, j, beta)) {
        fprintf(stderr, "error writing to file\n");
        exit(EXIT_FAILURE);
    }

#pragma omp threadprivate(xorshift_state)
#pragma omp parallel
    {
        for (int i = 0; i < omp_get_thread_num(); i++)
            jump();
    }

#pragma omp parallel for
    for (unsigned long i = 0; i < count; i++) {
        state_t lattice[TIME_LEN * SPINS_PER_STATE_T] = { 0 };
        initLattice(lattice);
        metropolis(lattice, hamiltonian(lattice, j, h_mu), j, h_mu, beta, iterations);
        // since writeState() only makes one call to fwrite, it should be thread-safe
        // this is actually only true on POSIX systems, linux and windows, but that is
        // basically all of the targets for this program
        if (writeState(data_file, lattice)) {
            fprintf(stderr, "error writing to file\n");
            exit(EXIT_FAILURE);
        }
        //puts("");
        //printLattice(lattice);
        //metropolis(lattice, hamiltonian(lattice, j, 0), j, 0, beta, iterations);
        //puts("zero applied field");
        //printLattice(lattice);

    }
    fclose(data_file);
}
