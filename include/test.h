#ifndef _TEST
#define _TEST

#include <openssl/rand.h>
#include <math.h>

void test_pcg();
double bench_pcg(size_t n, size_t c, size_t t);

void sample_a_and_a2(uint8_t *fft_a, uint32_t *fft_a2, size_t poly_size, size_t c);

#endif
