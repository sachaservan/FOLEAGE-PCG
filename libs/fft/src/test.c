#include <openssl/rand.h>
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>

#include <stdlib.h>
#include <time.h>

#include "fft.h"
#include "utils.h"

#define NUMVARS 16

double testFFT_uint64()
{
    size_t num_vars = NUMVARS;
    size_t num_coeffs = ipow(3, num_vars);

    uint64_t *coeffs = malloc(sizeof(uint64_t) * num_coeffs);
    RAND_bytes((uint8_t *)coeffs, sizeof(uint64_t) * num_coeffs);

    //************************************************
    printf("Benchmarking FFT evaluation with uint64_t packing \n");
    //************************************************

    clock_t t;
    t = clock();
    fft_recursive_uint64(coeffs, num_vars, num_coeffs / 3);
    t = clock() - t;
    double time_taken = ((double)t) / (CLOCKS_PER_SEC / 1000.0); // ms

    printf("FFT (uint64) eval time (total) %f ms\n", time_taken);

    free(coeffs);

    return time_taken;
}

double testFFT_uint32()
{
    size_t num_vars = NUMVARS;
    size_t num_coeffs = ipow(3, num_vars);

    uint32_t *coeffs = malloc(sizeof(uint32_t) * num_coeffs);
    RAND_bytes((uint8_t *)coeffs, sizeof(uint32_t) * num_coeffs);

    //************************************************
    printf("Benchmarking FFT evaluation with uint32_t packing \n");
    //************************************************

    clock_t t;
    t = clock();
    fft_recursive_uint32(coeffs, num_vars, num_coeffs / 3);
    t = clock() - t;
    double time_taken = ((double)t) / (CLOCKS_PER_SEC / 1000.0); // ms

    printf("FFT (uint32) eval time (total) %f ms\n", time_taken);

    free(coeffs);

    return time_taken;
}

double testFFT_uint8()
{
    size_t num_vars = NUMVARS;
    size_t num_coeffs = ipow(3, num_vars);
    uint8_t *coeffs = malloc(sizeof(uint8_t) * num_coeffs);
    RAND_bytes((uint8_t *)coeffs, sizeof(uint8_t) * num_coeffs);

    //************************************************
    printf("Benchmarking FFT evaluation without packing \n");
    //************************************************

    clock_t t;
    t = clock();
    fft_recursive_uint8(coeffs, num_vars, num_coeffs / 3);
    t = clock() - t;
    double time_taken = ((double)t) / (CLOCKS_PER_SEC / 1000.0); // ms

    printf("FFT (uint8) eval time (total) %f ms\n", time_taken);

    free(coeffs);

    return time_taken;
}

int main(int argc, char **argv)
{
    double time = 0;
    int testTrials = 5;

    printf("******************************************\n");
    printf("Testing FFT (uint8 packing)\n");
    for (int i = 0; i < testTrials; i++)
    {
        time += testFFT_uint8();
        printf("Done with trial %i of %i\n", i + 1, testTrials);
    }
    printf("******************************************\n");
    printf("DONE\n");
    printf("Avg time: %0.2f\n", time / testTrials);
    printf("******************************************\n\n");

    printf("******************************************\n");
    printf("Testing FFT (uint32 packing) \n");
    time = 0;
    for (int i = 0; i < testTrials; i++)
    {
        time += testFFT_uint32();
        printf("Done with trial %i of %i\n", i + 1, testTrials);
    }
    printf("******************************************\n");
    printf("DONE\n");
    printf("Avg time: %0.2f\n", time / testTrials);
    printf("******************************************\n\n");

    printf("******************************************\n");
    printf("Testing FFT (uint64 packing) \n");
    time = 0;
    for (int i = 0; i < testTrials; i++)
    {
        time += testFFT_uint64();
        printf("Done with trial %i of %i\n", i + 1, testTrials);
    }
    printf("******************************************\n");
    printf("DONE\n");
    printf("Avg time: %0.2f\n", time / testTrials);
    printf("******************************************\n\n");
}