#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
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

// apply fourier transform across the space dimension
// output should have TIME_LEN complex doubles allocated
void fourierTransformSpace(state_t *lattice, complex double *output, double p_n)
{
    for (int t = 0; t < TIME_LEN; t++) {
        complex double sum = CMPLX(0, 0);
        for (int x = 0; x < SPACE_LEN; x++)
            sum += CMPLX(getSpinAt(lattice, x, t), 0) * cexp(CMPLX(0, p_n * x));
        output[t] = sum;
    }
}

int main(int argc, char **argv)
{
    state_t lattice[TIME_LEN * SPACE_STATE_COUNT];
    int correlation_array[TIME_LEN * SPACE_LEN] = { 0 };

    if (argc != 3) {
        fprintf(stderr, "usage: %s <infile> <outfile>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    FILE *data_file = fopen(argv[1], "r");
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

    npy_array_t correlation_out = createNpyArrayNd('c', sizeof(complex double), 1, TIME_LEN);
    complex double *correlations = (complex double *)correlation_out.data;
    npy_array_t stddev_out = createNpyArrayNd('c', sizeof(complex double), 1, TIME_LEN);
    complex double *stddev = (complex double *)stddev_out.data;

    int state_counter;
    int n = 0;
    for (state_counter = 0; readState(data_file, lattice) == READ_SUCCESS; state_counter++) {
        complex double output[TIME_LEN];
        fourierTransformSpace(lattice, output, 2 * M_PI * n / SPACE_LEN);
        for (int i = 0; i < TIME_LEN; i++) {
            complex double temp = conj(output[0]) * output[i];
            correlations[i] += temp;
            // for now just use stddev to hold the sum of the squares,
            // we will need this later for variance calculation
            stddev[i] += temp * temp;
        }
    }

    for (int i = 0; i < TIME_LEN; i++) {
        correlations[i] /= state_counter;
        complex double mean = correlations[i];
        stddev[i] = csqrt((stddev[i] - state_counter * mean * mean) / state_counter) / sqrt(state_counter);
    }

    npy_array_list_t *array_head = npy_array_list_prepend(NULL, &stddev_out, "error");
    if (!array_head) {
        fprintf(stderr, "npy_array_list error\n");
        exit(EXIT_FAILURE);
    }
    array_head = npy_array_list_prepend(array_head, &correlation_out, "correlations");
    if (!array_head) {
        fprintf(stderr, "npy_array_list error\n");
        exit(EXIT_FAILURE);
    }
    
    if (npy_array_list_save(argv[2], array_head) != 2) {
        fprintf(stderr, "error saving array list\n");
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;

    //npy_array_save(argv[2], &correlation_out);

    // int state_counter;
    // for (state_counter = 0; readState(data_file, lattice) == READ_SUCCESS; state_counter++) {
    //     multiplyLatticeBy(lattice, getSpinAt(lattice, SPACE_LEN / 2, TIME_LEN / 2));
    //     addToCorrelationArray(lattice, correlation_array);
    // }

    // npy_array_t output = createNpyDoubleArrayNd(2, TIME_LEN, SPACE_LEN);

    // for (int t = 0; t < TIME_LEN; t++)
    //     for (int x = 0; x < SPACE_LEN; x++)
    //         ((double *)output.data)[t * SPACE_LEN + x] = (double)correlation_array[t * SPACE_LEN + x] / state_counter;

    // npy_array_save(argv[2], &output);
}
