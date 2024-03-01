#include <stdint.h>
#include <inttypes.h>
#include <immintrin.h>

#ifndef NUM
#define N (128U)
#else
#define N (NUM)
#endif

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

typedef __m256i mat_width_t;
#define BITS_PER_WIDTH ((uint32_t)sizeof(mat_width_t) * CHAR_BIT)
#define ROWS_PER_WIDTH (BITS_PER_WIDTH / N)

typedef mat_width_t sqr_mat[N / (ROWS_PER_WIDTH)];

extern sqr_mat R_i_lut[N];
extern sqr_mat L_i_lut[N];
extern size_t num_facts_of_P;
extern uint32_t factors_of_P[];

#if N == 32U
typedef uint32_t mat_row_t;
#define P ((uint64_t)UINT32_MAX)
#elif N == 64U
typedef uint64_t mat_row_t;
#define P ((uint64_t)UINT64_MAX)
#elif N == 128U
typedef __m128i mat_row_t;
#define P 0ULL

#else // default
typedef uint32_t mat_row_t;
#ifndef N
#define N (32U)
#endif
#define P ((uint64_t)UINT32_MAX)
#error "N must be (currently) either 32, 64 or 128"
#endif

#define MAT_ROWS_WIDTH(mat, i) ((mat)[(i)])

#define I (R_i_lut[0])

#endif // __CONSTANTS_H__