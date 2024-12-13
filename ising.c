#include "ising.h"

#include <errno.h>
#include <stdio.h>
#include <math.h>
#include "endian.h"
#include "popcount.h"
#include "npy_array/npy_array.h"

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

double hamiltonian(state_t *lattice, double j, double h_mu)
{
    // multiply horizontally
    int horizontal_energy = 0;
    int total_spins = 0;
    for (int i = 0; i < TIME_LEN; i++) {
        state_t last_carry = lattice[(i + 1) * SPACE_STATE_COUNT - 1] << (SPINS_PER_STATE_T - 1);
        for (int j = 0; j < SPACE_STATE_COUNT; j++) {
            state_t this_state = lattice[i * SPACE_STATE_COUNT + j];
            horizontal_energy += popcount(~((last_carry | this_state >> 1) ^ this_state));
            total_spins += popcount(this_state);
            last_carry = this_state << (SPINS_PER_STATE_T - 1);
        }
    }
    // Now, we need to subtract some constant from it in order to account for the fact that a 0 -> -1 and 1 -> 1
    horizontal_energy = 2 * horizontal_energy - TIME_LEN * SPACE_LEN;

    // multiply vertically
    int vertical_energy = 0;
    state_t *above_row = lattice + (TIME_LEN - 1) * SPACE_STATE_COUNT;
    for (int i = 0; i < TIME_LEN; i++) {
        state_t *current_row = lattice + i * SPACE_STATE_COUNT;
        for (int j = 0; j < SPACE_STATE_COUNT; j++)
            vertical_energy += popcount(~(current_row[j] ^ above_row[j]));
        above_row = current_row;
    }
    // Same deal here
    vertical_energy = 2 * vertical_energy - TIME_LEN * SPACE_LEN;

    return -j * (horizontal_energy + vertical_energy) - h_mu * (2 * total_spins - TIME_LEN * SPACE_LEN);
}

double hamiltonianDebug(state_t *lattice, double j, double h_mu)
{
    int total_energy = 0;
    for (int x = 0; x < SPACE_LEN; x++) {
        for (int t = 0; t < TIME_LEN; t++) {
            spin_t this_spin = getSpinAt(lattice, x, t);
            total_energy -= j * (this_spin * getSpinAt(lattice, x - 1, t) + this_spin * getSpinAt(lattice, x, t - 1))
                + h_mu * this_spin;
        }
    }
    return total_energy;
}

double calculateEnergyChange(state_t *lattice, double j, double h_mu, int x, int t)
{
    int up     = getSpinAt(lattice, x, t - 1);
    int down   = getSpinAt(lattice, x, t + 1);
    int left   = getSpinAt(lattice, x - 1, t);
    int right  = getSpinAt(lattice, x + 1, t);
    int center = getSpinAt(lattice, x, t);
    
    // first, calculate its current contribution
    int current = center * (up + down + left + right);
    // the new contribution is just -current
    return 2.0 * (j * (double)current + h_mu * center);
}

double metropolis(state_t *lattice, double energy, double j, double h_mu, double beta, int iterations)
{
    for (int i = 0; i < iterations; i++) {
        // first, pick a random point in spacetime
        int x = randomInt(0, SPACE_LEN);
        int t = randomInt(0, TIME_LEN);

        double delta = calculateEnergyChange(lattice, j, h_mu, x, t);

        if (uniformFloat() <= exp(-beta * delta)) {
            flipSpinAt(lattice, x, t);
            energy += delta;
        }
    }
    return energy;
}
