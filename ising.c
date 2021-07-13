#include <limits.h>
#include <stdint.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>
#include "npy_array/npy_array.h"

typedef uint64_t state_t; // several spins in a compressed format for fast batch operations
typedef int spin_t;       // the type used to represent a single spin

#define SPINS_PER_STATE_T (sizeof(state_t) * CHAR_BIT)                                // the number of spins stored in one state_t variable, aka the number of bits in state_t
#define TIME_LEN 64                                                                   // size of lattice in the time dimension
#define SPACE_LEN 64                                                                  // size of lattice in the space dimension
#define SPACE_STATE_COUNT ((SPACE_LEN + SPINS_PER_STATE_T - 1) / SPINS_PER_STATE_T)   // number of state_t elements in the space dimension
#define SPACE_REMAINDER ((SPINS_PER_STATE_T * SPACE_STATE_COUNT - 1) % SPACE_LEN + 1)
//       space
//      *------>
// time | 0, 1
//      | 2, 3
//      v

// just some random bytes I grabbed off RANDOM.org
uint64_t xorshift_state[4] = { 
    0x5959c5803c47d89b, 0x21eeb5e4b0c8ae73,
    0xdb8c526a43d863fe, 0x1cae287ffe7ad6fd 
};

uint64_t bitRoll64(uint64_t a, int k)
{
    return (a << k) | (a >> (64 - k));
}

// I will use this random number generator (xorshiro256**) to fill the lattice
uint64_t xorshift256()
{
    uint64_t result = bitRoll64(xorshift_state[1] * 5, 7) * 9;
    uint64_t temp = xorshift_state[1] << 17;

    xorshift_state[2] ^= xorshift_state[0];
    xorshift_state[3] ^= xorshift_state[1];
    xorshift_state[1] ^= xorshift_state[2];
    xorshift_state[0] ^= xorshift_state[3];

    xorshift_state[2] ^= temp;
    xorshift_state[3] ^= bitRoll64(xorshift_state[3], 45);
    
    return result;
}

int randomInt(int lower, int upper)
{
    int modulo = upper - lower;
    uint64_t random = xorshift256();
    uint64_t unbiased_max = 0xffffffffffffffff - (0xffffffffffffffff % modulo + 1) % modulo;
    while (__builtin_expect(random > unbiased_max, 0))
        random = xorshift256();
    return lower + random % modulo;
}

double uniformFloat()
{
    return (xorshift256() >> 11) * 0x1.0p-53;
}

// hot start the lattice
void initLattice(state_t *lattice)
{
    for (int i = 0; i < TIME_LEN * SPACE_STATE_COUNT; i++) {
        lattice[i] = (state_t)xorshift256();
    }
}

// takes the lattice, a position x, and a time t, and returns the spin at that point as either +1 or -1
spin_t getSpinAt(state_t *lattice, int x, int t)
{
    // this could be optimized significantly if we could assume that SPACE_LEN and TIME_LEN are powers of 2
    // not sure if that is ok to assume though 
    x = x % SPACE_LEN + SPACE_LEN * (x < 0); // boundary condition is periodic so we need to take the modulus
    t = t % TIME_LEN + TIME_LEN * (t < 0);
    int bit = (lattice[x / SPINS_PER_STATE_T + t * SPACE_STATE_COUNT] >> (x % SPINS_PER_STATE_T)) & 1;
    return 2 * bit - 1;
}

// takes the lattice, a position x, and a time t, and returns the spin at that point as either +1 or -1
void setSpinAt(state_t *lattice, int x, int t, int spin)
{
    // this could be optimized significantly if we could assume that SPACE_LEN and TIME_LEN are powers of 2
    // not sure if that is ok to assume though 
    x = x % SPACE_LEN + SPACE_LEN * (x < 0); // boundary condition is periodic so we need to take the modulus
    t = t % TIME_LEN + TIME_LEN * (t < 0);
    state_t bit = 1 << (x % SPINS_PER_STATE_T);
    state_t value = bit * ((spin + 1) / 2);
    lattice[x / SPINS_PER_STATE_T + t * SPACE_STATE_COUNT] &= ~bit;
    lattice[x / SPINS_PER_STATE_T + t * SPACE_STATE_COUNT] |= value;
}

void flipSpinAt(state_t *lattice, int x, int t)
{
    // this could be optimized significantly if we could assume that SPACE_LEN and TIME_LEN are powers of 2
    // not sure if that is ok to assume though 
    x = x % SPACE_LEN + SPACE_LEN * (x < 0); // boundary condition is periodic so we need to take the modulus
    t = t % TIME_LEN + TIME_LEN * (t < 0);
    lattice[x / SPINS_PER_STATE_T + t * SPACE_STATE_COUNT] ^= (state_t)1 << (x % SPINS_PER_STATE_T);
}

void printLattice(state_t *lattice)
{
    for (int i = 0; i < TIME_LEN; i++) {
        printf("< ");
        for (int j = 0; j < SPACE_LEN - 1; j++) {
            // write either +1 or -1 depending on the spin
            // ASCII 43: '+' ASCII 44: ',' ASCII 45: '-'
            printf("%c ", ',' - getSpinAt(lattice, j, i));
        }
        printf("%c >\n", ',' - getSpinAt(lattice, SPACE_LEN - 1, i));
    }
}

double hamiltonian(state_t *lattice, double j)
{
    //printf("%d\n", SPACE_REMAINDER);
    // multiply horizontally
    int horizontal_energy = 0;
    for (int i = 0; i < TIME_LEN; i++) {
        state_t last_carry = lattice[(i + 1) * SPACE_STATE_COUNT - 1] << (SPINS_PER_STATE_T - 1);
        //printf("%lx\n", last_carry);
        for (int j = 0; j < SPACE_STATE_COUNT - 1; j++) {
            state_t this_state = lattice[i * SPACE_STATE_COUNT + j];
            horizontal_energy += __builtin_popcountl(~((last_carry | this_state >> 1) ^ this_state));
            last_carry = this_state << (SPINS_PER_STATE_T - 1);
        }
        state_t this_state = lattice[(i + 1) * SPACE_STATE_COUNT - 1];
        //printf("this_state: %lx\n", last_carry >> (SPINS_PER_STATE_T - SPACE_REMAINDER) | this_state >> 1);
        horizontal_energy += __builtin_popcountl(~((last_carry >> (SPINS_PER_STATE_T - SPACE_REMAINDER) | this_state >> 1) ^ this_state) & (state_t)-1 >> (SPINS_PER_STATE_T - SPACE_REMAINDER));
    }
    // Now, we need to subtract some constant from it in order to account for the fact that a 0 -> -1 and 1 -> 1
    horizontal_energy = 2 * horizontal_energy - TIME_LEN * SPACE_LEN;
    
    // multiply vertically
    int vertical_energy = 0;
    state_t *above_row = lattice + (TIME_LEN - 1) * SPACE_STATE_COUNT;
    for (int i = 0; i < TIME_LEN; i++) {
        state_t *current_row = lattice + i * SPACE_STATE_COUNT;
        for (int j = 0; j < SPACE_STATE_COUNT - 1; j++)
            vertical_energy += __builtin_popcountl(~(current_row[j] ^ above_row[j]));
        vertical_energy += __builtin_popcountl(~(current_row[SPACE_STATE_COUNT - 1] ^ above_row[SPACE_STATE_COUNT - 1]) & ((state_t)-1 >> (SPINS_PER_STATE_T - SPACE_REMAINDER)));
        above_row = current_row;
    }
    // Same deal here
    vertical_energy = 2 * vertical_energy - TIME_LEN * SPACE_LEN;

    return -j * (horizontal_energy + vertical_energy);
}

double hamiltonianDebug(state_t *lattice, double j)
{
    int total_energy = 0;
    for (int x = 0; x < SPACE_LEN; x++) {
        for (int t = 0; t < TIME_LEN; t++) {
            spin_t this_spin = getSpinAt(lattice, x, t);
            total_energy += this_spin * getSpinAt(lattice, x - 1, t) + this_spin * getSpinAt(lattice, x, t - 1);
        }
    }
    return -j * total_energy;
}

double calculateEnergyChange(state_t *lattice, double j, int x, int t)
{
    int up     = getSpinAt(lattice, x, t - 1);
    int down   = getSpinAt(lattice, x, t + 1);
    int left   = getSpinAt(lattice, x - 1, t);
    int right  = getSpinAt(lattice, x + 1, t);
    int center = getSpinAt(lattice, x, t);
    
    // first, calculate its current contribution
    int current = center * (up + down + left + right);
    // the new contribution is just -current
    return 2.0 * j * (double)current;
}

double metropolis(state_t *lattice, double energy, double j, double beta, int iterations)
{
    for (int i = 0; i < iterations; i++) {
        // first, pick a random point in spacetime
        int x = randomInt(0, SPACE_LEN);
        int t = randomInt(0, TIME_LEN);

        double delta = calculateEnergyChange(lattice, j, x, t);

        if (uniformFloat() <= exp(-beta * delta)) {
            flipSpinAt(lattice, x, t);
            energy += delta;
        }
    }
    return energy;
}

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
