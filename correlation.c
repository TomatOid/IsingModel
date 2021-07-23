#include <stdio.h>
#include <stdlib.h>
#include "record.h"
#include "ising.h"

void multiplyLatticeBy(state_t *lattice, int n)
{
    state_t xor_by = (state_t)-1;
    if (n == 1)
        return;
    
    for (int t = 0; t < TIME_LEN; t++)
        for (int x = 0; x < SPACE_STATE_COUNT; x++)
            lattice[t * SPACE_STATE_COUNT + x] ^= xor_by;
}

void addToCorrelationArray(state_t *lattice, int *correlation_array)
{
    for (int t = 0; t < TIME_LEN; t++) {
        for (int x = 0; x < SPACE_STATE_COUNT; x++) {
            state_t spin_set = lattice[t * SPACE_STATE_COUNT + x];

            for (int bit = 0; bit < SPINS_PER_STATE_T; bit++, spin_set >>= 1)
                correlation_array[t * SPACE_LEN + x * SPINS_PER_STATE_T + bit] +=
                    2 * (int)(spin_set & 1) - 1;
        }
    }
}

int main()
{
    state_t lattice[TIME_LEN * SPACE_STATE_COUNT];
    int correlation_array[TIME_LEN * SPACE_LEN] = { 0 };

    FILE *data_file = fopen("states.dat", "r");
    if (!data_file) {
        fputs("error opening data file", stderr);
        exit(EXIT_FAILURE);
    }
    
    double j, beta;
    int error_code = readHeader(data_file, &j, &beta);
    if (error_code != READ_SUCCESS) {
        printf("error with header, code %d\n", error_code);
        exit(EXIT_FAILURE);
    }

    printf("j: %f, beta: %f\n", j, beta);

    int state_counter;
    for (state_counter = 0; readState(data_file, lattice) == READ_SUCCESS; state_counter++) {
        multiplyLatticeBy(lattice, getSpinAt(lattice, 0, 0));
        addToCorrelationArray(lattice, correlation_array);
    }

    for (int t = 0; t < TIME_LEN; t++) {
        for (int x = 0; x < SPACE_LEN - 1; x++)
            printf("%f, ", ((double)correlation_array[t * SPACE_LEN + x]) / state_counter);
        printf("%f\n", ((double)correlation_array[(t + 1) * SPACE_LEN - 1]) / state_counter);
    }
}
