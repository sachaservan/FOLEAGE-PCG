#ifndef _FFT
#define _FFT

#include <string.h>
#include <stdint.h>

typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

// FFT for (up to) 32 polynomials over F4
void fft_recursive_uint64(
    uint64_t *coeffs,
    const size_t num_vars,
    const size_t num_coeffs);

// FFT for (up to) 16 polynomials over F4
void fft_recursive_uint32(
    uint32_t *coeffs,
    const size_t num_vars,
    const size_t num_coeffs);

// FFT for (up to) 8 polynomials over F4
void fft_recursive_uint16(
    uint16_t *coeffs,
    const size_t num_vars,
    const size_t num_coeffs);

// FFT for (up to) 4 polynomials over F4
void fft_recursive_uint8(
    uint8_t *coeffs,
    const size_t num_vars,
    const size_t num_coeffs);

// Compute base^exp without the floating-point precision
// errors of the built-in pow function.
static inline size_t ipow(size_t base, size_t exp)
{
    if (exp == 1)
        return base;

    if (exp == 0)
        return 1;

    size_t result = 1;
    while (1)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

#endif
