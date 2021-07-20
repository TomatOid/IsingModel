#pragma once
#include <limits.h>
#include <stdint.h>

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

extern uint64_t xorshift_state[4];

// xorshiro256** PRNG
uint64_t xorshift256();

// generates a random number on the interval [lower, upper)
int randomInt(int lower, int upper);

// generates a double-precision floating point number with a
// uniform distribution from 0 to 1
double uniformFloat();

// uses the xorshiro256** PRNG to hot start the lattice
void initLattice(state_t *lattice);

// takes the lattice, a position x, and a time t, and returns
// the spin at that point as either +1 or -1
spin_t getSpinAt(state_t *lattice, int x, int t);

// flips a spin at some position x and time t
void flipSpinAt(state_t *lattice, int x, int t);

void printLattice(state_t *lattice);

double hamiltonian(state_t *lattice, double j);

double calculateEnergyChange(state_t *lattice, double j, int x, int t);

double metropolis(state_t *lattice, double energy, double j, double beta, int iterations);
