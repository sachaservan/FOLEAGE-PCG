#ifndef _F4OPS
#define _F4OPS

typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

// Samples a non-zero element of F4
static uint8_t rand_f4x()
{
    uint8_t t;
    unsigned char rand_byte;

    // loop until we have two bits where at least one is non-zero
    while (1)
    {
        RAND_bytes(&rand_byte, 1);
        t = 0;
        t |= rand_byte & 1;
        t = t << 1;
        t |= (rand_byte >> 1) & 1;

        if (t != 0 && t != 4)
            return t;
    }
}

// Multiplies two elements of F4 (optionally: 4 elements packed into uint8_t)
// and returns the result.
static uint8_t mult_f4(uint8_t a, uint8_t b)
{
    uint128_t tmp = ((a & 0b10) & (b & 0b10));
    uint8_t res = tmp ^ ((a & 0b10) & ((b & 0b01) << 1) ^ (((a & 0b01) << 1) & (b & 0b10)));
    res |= ((a & 0b01) & (b & 0b01)) ^ (tmp >> 1);
    return res;
}

// Multiplies two packed matrices of F4 elements column-by-column.
// Note that here the "columns" are packed into an element of uint8_t
// resulting in a matrix with 4 columns.
static void multiply_fft_8(
    uint8_t *a_poly,
    uint8_t *b_poly,
    uint8_t *res_poly,
    size_t poly_size)
{
    const uint8_t pattern = 0xaa;
    uint8_t mask_h = pattern;     // 0b101010101010101001010
    uint8_t mask_l = mask_h >> 1; // 0b010101010101010100101

    uint8_t tmp;
    uint8_t a_h, a_l, b_h, b_l;

    for (size_t i = 0; i < poly_size; i++)
    {
        // multiplication over F4
        a_h = (a_poly[i] & mask_h);
        a_l = (a_poly[i] & mask_l);
        b_h = (b_poly[i] & mask_h);
        b_l = (b_poly[i] & mask_l);

        tmp = (a_h & b_h);
        res_poly[i] = tmp ^ (a_h & (b_l << 1));
        res_poly[i] ^= ((a_l << 1) & b_h);
        res_poly[i] |= a_l & b_l ^ (tmp >> 1);
    }
}

// Multiplies two packed matrices of F4 elements column-by-column.
// Note that here the "columns" are packed into an element of uint16_t
// resulting in a matrix with 8 columns.
static void multiply_fft_16(
    uint16_t *a_poly,
    uint16_t *b_poly,
    uint16_t *res_poly,
    size_t poly_size)
{
    const uint16_t pattern = 0xaaaa;
    uint16_t mask_h = pattern;     // 0b101010101010101001010
    uint16_t mask_l = mask_h >> 1; // 0b010101010101010100101

    uint16_t tmp;
    uint16_t a_h, a_l, b_h, b_l;

    for (size_t i = 0; i < poly_size; i++)
    {
        // multiplication over F4
        a_h = (a_poly[i] & mask_h);
        a_l = (a_poly[i] & mask_l);
        b_h = (b_poly[i] & mask_h);
        b_l = (b_poly[i] & mask_l);

        tmp = (a_h & b_h);
        res_poly[i] = tmp ^ (a_h & (b_l << 1));
        res_poly[i] ^= ((a_l << 1) & b_h);
        res_poly[i] |= a_l & b_l ^ (tmp >> 1);
    }
}

// Multiplies two packed matrices of F4 elements column-by-column.
// Note that here the "columns" are packed into an element of uint32_t
// resulting in a matrix with 16 columns.
static void multiply_fft_32(
    uint32_t *a_poly,
    uint32_t *b_poly,
    uint32_t *res_poly,
    size_t poly_size)
{
    const uint32_t pattern = 0xaaaaaaaa;
    uint32_t mask_h = pattern;     // 0b101010101010101001010
    uint32_t mask_l = mask_h >> 1; // 0b010101010101010100101

    uint32_t tmp;
    uint32_t a_h, a_l, b_h, b_l;

    for (size_t i = 0; i < poly_size; i++)
    {
        // multiplication over F4
        a_h = (a_poly[i] & mask_h);
        a_l = (a_poly[i] & mask_l);
        b_h = (b_poly[i] & mask_h);
        b_l = (b_poly[i] & mask_l);

        tmp = (a_h & b_h);
        res_poly[i] = tmp ^ (a_h & (b_l << 1));
        res_poly[i] ^= ((a_l << 1) & b_h);
        res_poly[i] |= a_l & b_l ^ (tmp >> 1);
    }
}

// Multiplies two packed matrices of F4 elements column-by-column.
// Note that here the "columns" are packed into an element of uint64_t
// resulting in a matrix with 32 columns.
static void multiply_fft_64(
    uint64_t *a_poly,
    uint64_t *b_poly,
    uint64_t *res_poly,
    size_t poly_size)
{
    const uint64_t pattern = 0xaaaaaaaaaaaaaaaa;
    uint64_t mask_h = pattern;     // 0b101010101010101001010
    uint64_t mask_l = mask_h >> 1; // 0b010101010101010100101

    uint64_t tmp;
    uint64_t a_h, a_l, b_h, b_l;

    for (size_t i = 0; i < poly_size; i++)
    {
        // multiplication over F4
        a_h = (a_poly[i] & mask_h);
        a_l = (a_poly[i] & mask_l);
        b_h = (b_poly[i] & mask_h);
        b_l = (b_poly[i] & mask_l);

        tmp = (a_h & b_h);
        res_poly[i] = tmp ^ (a_h & (b_l << 1));
        res_poly[i] ^= ((a_l << 1) & b_h);
        res_poly[i] |= a_l & b_l ^ (tmp >> 1);
    }
}

#endif