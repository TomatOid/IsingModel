#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "record.h"
#include "ising.h"

// apply fourier transform across the space dimension
// output should have TIME_LEN complex doubles allocated
void fourierTransformSpace(state_t *lattice, complex double *output, unsigned k)
{
    double p_n = 2 * M_PI * (double)k / SPACE_LEN;

    static complex double dft_precomp_table[SPACE_LEN][2 * SPACE_LEN] = {{ 0 }};
    static int initialized[SPACE_LEN] = { 0 };
    complex double *dft_precomp = NULL;
    if (k < SPACE_LEN)
        dft_precomp = &dft_precomp_table[k][0];
    else
        return;

    if (!initialized[k]) {
        for (int x = 0; x < 2 * SPACE_LEN; x += 2) {
            complex double vec = cexp(CMPLX(0, p_n * (x / 2)));
            dft_precomp[x] = -vec;
            dft_precomp[x + 1] = vec;
        }
        initialized[k] = 1;
    }

    for (int t = 0; t < TIME_LEN; t++) {
        complex double sum = CMPLX(0, 0);
        for (int x = 0; x < SPACE_STATE_COUNT; x++) {
            state_t spin_set = lattice[t * SPACE_STATE_COUNT + x];
            
            for (int bit = 0; bit < SPINS_PER_STATE_T; bit++, spin_set >>= 1)
                sum += dft_precomp[2 * (x * SPINS_PER_STATE_T + bit) + (spin_set & 1)];
        }
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

    npy_array_t correlation_out = createNpyArrayNd('c', sizeof(complex double), 2, SPACE_LEN / 2, TIME_LEN);
    complex double *correlations = (complex double *)correlation_out.data;
    npy_array_t stddev_out = createNpyArrayNd('f', sizeof(double), 2, SPACE_LEN / 2, TIME_LEN);
    double *stddev = (double *)stddev_out.data;

    int state_counter;
    for (state_counter = 0; readState(data_file, lattice) == READ_SUCCESS; state_counter++) {
        for (int n = 0; n < SPACE_LEN / 2; n++) {
            complex double output[TIME_LEN];
            fourierTransformSpace(lattice, output, n);
            for (int i = 0; i < TIME_LEN; i++) {
                complex double temp = conj(output[0]) * output[i];

                // welford's algo
                complex double delta = temp - correlations[n * TIME_LEN + i];
                correlations[n * TIME_LEN + i] += delta / (state_counter + 1);
                complex double delta2 = temp - correlations[n * TIME_LEN + i];
                stddev[n * TIME_LEN + i] += cabs(delta) * cabs(delta2);
            }
        }
    }

    for (int n = 0; n < SPACE_LEN / 2; n++) {
        double norm = cabs(correlations[n * TIME_LEN]);
        for (int i = 0; i < TIME_LEN; i++) {
            correlations[n * TIME_LEN + i] /= norm;
            stddev[n * TIME_LEN + i] = sqrt(stddev[n * TIME_LEN + i] / state_counter) / norm / sqrt(state_counter);
        }
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
}
