#ifndef _UTILS
#define _UTILS

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