#ifdef __SQUARE_MATRIX_IMPLEMENTATION__

#ifndef __SQUARE_MATRIX_SIZE__
#error "[ERRROR] you must define a matrix size before including this header !!!"
#endif

#ifndef __SQUARE_MATRIX_NAMESPACE__
#error "[ERROR] you must define a namespace before including this header !!!"
#endif

#ifndef __SQUARE_MATRIX_TYPE__
#error "[ERROR] you must define a matrix type before including this header !!!"
#endif

#ifndef __SQUARE_MATRIX_ROW_TYPE__
#error "[ERROR] you must define a matrix row type before including this header !!!"
#endif

#ifndef __SQUARE_MATRIX_ROWS_PER_WIDTH__
#if __SQUARE_MATRIX_ENABLE_WARNINGS__
#warning "[WARNING] you should define the number of rows per width before including this header !!!"
#endif
#define __SQUARE_MATRIX_ROWS_PER_WIDTH__ ROWS_PER_WIDTH /*default*/
#endif

#ifndef __SQUARE_MATRIX_GET_MAT__
#if __SQUARE_MATRIX_ENABLE_WARNINGS__
#warning "[WARNING] you should define a way to get a matrix before including this header !!!"
#endif
#define __SQUARE_MATRIX_GET_MAT__ GET_MAT
#endif

#define XCONC(A, B) A##B
#define CONC(A, B) XCONC(A, B)

#define FN_NAME(str) CONC(__SQUARE_MATRIX_NAMESPACE__, str)
#define SUBMAT_FN_NAME(str) CONC(__SQUARE_MATRIX_SUBMAT_NAMESPACE__, str)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define _mm256_set_zero() _mm256_set1_epi64x(0)
#if __SQUARE_MATRIX_SIZE__ == 32U
#define _mm256_set1(uint32_t) _mm256_set1_epi32(uint32_t)
#define _mm256_maskload(const_int_ptr, __m256i) _mm256_maskload_epi32((const int *)(const_int_ptr), __m256i)
#define _mm256_slli(__m256i, int) _mm256_slli_epi32(__m256i, int)
#define SHIFTS_OFF_VALS 0, 1, 2, 3, 4, 5, 6, 7
#elif __SQUARE_MATRIX_SIZE__ == 64U
#define _mm256_set1(uint64_t) _mm256_set1_epi64x(uint64_t)
#define _mm256_setr(...) _mm256_setr_epi64x(__VA_ARGS__)
#define _mm256_maskload(const_int_ptr, __m256i) _mm256_maskload_epi64((const long long int *)(const_int_ptr), __m256i)
#define _mm256_slli(__m256i, int) _mm256_slli_epi64(__m256i, int)
#define SHIFTS_OFF_VALS 0, 1, 2, 3
#elif __SQUARE_MATRIX_SIZE__ == 128U
#define _mm256_set1(__m128i) _mm256_broadcastsi128_si256(__m128i)
#define _mm256_slli(__m256i, int) (_mm256_permutevar8x32_epi32(_mm256_slli_epi32(__m256i, int % 32), _mm256_add_epi32(_mm256_setr_epi32(0, 0, 0, 0, 4, 4, 4, 4), _mm256_set1_epi32(int / 32))))
#define _mm256_maskload(const_int_ptr, __m256i_mask) _mm256_maskload_epi32((const int *)(const_int_ptr), __m256i_mask)
#else
#error "N must (currently) either 32, 64 or 128"
#endif

int FN_NAME(add)(__SQUARE_MATRIX_TYPE__ *res, __SQUARE_MATRIX_TYPE__ a, __SQUARE_MATRIX_TYPE__ b)
{
#ifndef __SQUARE_MATRIX_MULT_SUBMAT__

  for (uint8_t i = 0; i < (__SQUARE_MATRIX_SIZE__ / __SQUARE_MATRIX_ROWS_PER_WIDTH__); ++i)
    MAT_ROWS_WIDTH(*res, i) = _mm256_xor_si256(MAT_ROWS_WIDTH(a, i), MAT_ROWS_WIDTH(b, i));
  return 0;

#else

  for (uint8_t i = 0; i < __SQUARE_MATRIX_SUBROW_COUNT__; ++i)
    SUBMAT_FN_NAME(add)
    ()

#endif
}

int FN_NAME(mult)(__SQUARE_MATRIX_TYPE__ *restrict res, __SQUARE_MATRIX_TYPE__ a, __SQUARE_MATRIX_TYPE__ b)
{
#if __SQUARE_MATRIX_SIZE__ <= 64U || !defined(__SQUARE_MATRIX_MULT_SUBMAT__)

#if __SQUARE_MATRIX_SIZE__ == 128U
  __m256i mask = _mm256_set1_epi32(0x80000000U);
#else
  __m256i mask = _mm256_set1(((__SQUARE_MATRIX_ROW_TYPE__)1) << (__SQUARE_MATRIX_SIZE__ - 1)); // 0b100000...000
#endif

  for (uint32_t i = 0; i < (__SQUARE_MATRIX_SIZE__ / __SQUARE_MATRIX_ROWS_PER_WIDTH__); ++i)
  {
    MAT_ROWS_WIDTH(*res, i) = _mm256_set_zero();

    for (uint32_t j = 0; j < __SQUARE_MATRIX_SIZE__; ++j)
    {
      __m256i shifted_a = _mm256_slli(MAT_ROWS_WIDTH(a, i), j);

      __m256i masked_a = _mm256_and_si256(shifted_a, mask); // the highest bit tells if the row from b should be xored with res

      __m256i v_b_j = _mm256_set1(((__SQUARE_MATRIX_ROW_TYPE__ *)b)[j]);

      __m256i vals = _mm256_maskload((mat_row_t *)&v_b_j, masked_a);
      MAT_ROWS_WIDTH(*res, i) = _mm256_xor_si256(MAT_ROWS_WIDTH(*res, i), vals);
    }
  }

#elif defined(__SQUARE_MATRIX_MULT_STRASSEN__)

#else

#endif

  return 0;
}

int FN_NAME(pow_2_to_n)(__SQUARE_MATRIX_TYPE__ *restrict res, __SQUARE_MATRIX_TYPE__ a, uint32_t n)
{
  __SQUARE_MATRIX_TYPE__ tmp = {0};
  memcpy(&tmp, a, sizeof(tmp));
  for (uint32_t i = 0; i < n; ++i)
  {
    FN_NAME(mult)
    (res, tmp, tmp);
    memcpy(&tmp, res, sizeof(tmp));
  }
  return 0;
}

int FN_NAME(bin_exp)(__SQUARE_MATRIX_TYPE__ *restrict res, __SQUARE_MATRIX_TYPE__ a, size_t n)
{
  memcpy(res, &I, sizeof(*res));

  __SQUARE_MATRIX_TYPE__ a_ = {0};
  memcpy(&a_, a, sizeof(a_));
  while (n > 0)
  {
    if (n & 1)
    {
      // res = res * a
      __SQUARE_MATRIX_TYPE__ tmp = {0};
      FN_NAME(mult)
      (&tmp, *res, a_);
      memcpy(res, &tmp, sizeof(*res));
    }
    // a = a * a;
    __SQUARE_MATRIX_TYPE__ tmp = {0};
    FN_NAME(mult)
    (&tmp, a_, a_);
    memcpy(&a_, &tmp, sizeof(a_));

    n >>= 1;
  }

  return 0;
}

int FN_NAME(print)(__SQUARE_MATRIX_TYPE__ val)
{
  for (uint8_t j = 0; j < __SQUARE_MATRIX_SIZE__; ++j)
  {
    for (uint8_t i = 0; i < __SQUARE_MATRIX_SIZE__; ++i)
      printf("%c ", __SQUARE_MATRIX_GET_MAT__(val, i, j) ? '#' : (i % 32 == 0 && j % 32 == 0 ? 'o' : '.'));
    // putchar('\t');
    // for (uint8_t i = 0; i < SUBROWS_PER_ROW; ++i)
    //   printf("%08" PRIx32 " ", MAT_SUBROW(val, SUBROWS_PER_ROW * j + i));
    putchar('\n');
  }

  return 0;
}

#undef _mm256_set_zero
#undef _mm256_set1
#undef _mm256_maskload
#undef _mm256_slli
#undef SHIFTS_OFF_VALS

#undef __SQUARE_MATRIX_IMPLEMENTATION__
#undef __SQUARE_MATRIX_SIZE__
#undef __SQUARE_MATRIX_TYPE__
#undef __SQUARE_MATRIX_ROW_TYPE__
#undef __SQUARE_MATRIX_NAMESPACE__
#undef __SQUARE_MATRIX_ROWS_PER_WIDTH__
#undef __SQUARE_MATRIX_GET_MAT__

#endif // __SQUARE_MATRIX_IMPLEMENTATION__