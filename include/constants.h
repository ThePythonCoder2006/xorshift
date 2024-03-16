#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <immintrin.h>

#ifndef NUM
#define N (128U)
#else
#define N (NUM)
#endif

typedef __m256i mat_width_t;
#define BITS_PER_WIDTH ((uint32_t)sizeof(mat_width_t) * CHAR_BIT)
#define ROWS_PER_WIDTH (BITS_PER_WIDTH / N)

typedef mat_width_t sqr_mat[N / (ROWS_PER_WIDTH)];

#if N == 32U
typedef uint32_t mat_row_t;
#define P ((uint64_t)UINT32_MAX)
#define LUT_FNAME "data/lut_32.data"
#elif N == 64U
typedef uint64_t mat_row_t;
#define P ((uint64_t)UINT64_MAX)
#define LUT_FNAME "data/lut_64.data"
#elif N == 128U
typedef __m128i mat_row_t;
#define P 0ULL
#define LUT_FNAME "data/lut_32.data"
#else // default
typedef uint32_t mat_row_t;
#ifndef N
#define N (32U)
#endif
#define P ((uint64_t)UINT32_MAX)
#error "N must be (currently) either 32, 64 or 128"
#endif

#define MAT_ROWS_WIDTH(mat, i) ((mat)[(i)])
#define MAT_ROW(mat, i) ((mat_row_t *)&MAT_ROWS_WIDTH(mat, i / ROWS_PER_WIDTH))[i % ROWS_PER_WIDTH]

#define CLEAR_BIT(dst, i) ((dst) & (~(0x1 << (i))))

#define SET_BIT_TO_VAL(dst, i, val) \
  do                                \
  {                                 \
    (dst) = CLEAR_BIT((dst), (i));  \
    (dst) |= (val) << (i);          \
  } while (0)
// val MUST only have only its lowest bit set
#define SET_MAT(mat, i, j, val) SET_BIT_TO_VAL(MAT_ROWS_WIDTH(mat, j), N - (i)-1, val)

#define GET_BIT(dst, i) (((dst) >> (i)) & 0x1)
#if N == 128U
#define MAT_SUBROW(mat, i) ((subrow_t *)&(MAT_ROW((mat), (i) / SUBROWS_PER_ROW)))[(i) % SUBROWS_PER_ROW]
#define GET_MAT(mat, i, j) GET_BIT(MAT_SUBROW(mat, j *SUBROWS_PER_ROW + i / 32), 32 - (i)-1)
#else
#define GET_MAT(mat, i, j) GET_BIT(MAT_ROW(mat, j), N - (i)-1)
#endif // N == 128U

#if N > 64U

#ifndef SNUM
#define SUBN (32U)
#else
#define SUBN (SNUM)
#endif // SNUM

#if SUBN == 32U
typedef uint32_t subrow_t;
#elif SUBN == 64U
typedef uint64_t subrow_t;
#else // default
#error "only 32 and 64 bit subrow sizes are supported !!"
#endif // SUBN == 32U

#define SUBROWS_PER_WIDTH (sizeof(mat_width_t) / sizeof(subrow_t))
#define SUBROWS_PER_ROW (sizeof(mat_row_t) / sizeof(subrow_t))
typedef __m256i subwidth_t;

#define SUBROWS_PER_SUBWIDTH (sizeof(subwidth_t) / sizeof(subrow_t))
typedef subwidth_t submat[SUBN / (SUBROWS_PER_SUBWIDTH)];

#define SUBMAT_ROW(submat, i) (((subrow_t *)&MAT_ROWS_WIDTH(submat, i / SUBROWS_PER_SUBWIDTH))[i % SUBROWS_PER_SUBWIDTH])

#define GET_SUBMAT(mat, i, j) GET_BIT(MAT_SUBROW(mat, j), N - (i)-1)

#else // N > 64
#define SUBN (N)
typedef sqr_mat submat;
#endif // N > 64

extern submat R_i_lut[N];
extern submat L_i_lut[N];

#define I (R_i_lut[0])
