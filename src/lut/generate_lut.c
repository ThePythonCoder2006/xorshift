#define __USE_MINGW_ANSI_STDIO 1

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "constants.h"

#undef I

#if N > 64
#undef N
#define N (SUBN)
#endif

int sqr_mat_R_shift_n(submat *ret, uint8_t n);
int sqr_mat_L_shift_n(submat *ret, uint8_t n);
int sqr_mat_I(submat *ret);
int write_sqr_mat(FILE *f, submat val);

int main(int argc, char **argv)
{
  --argc, ++argv; // remove program name

  // default name
  char fname_default[] = LUT_FNAME;
  char *fname = fname_default;
  if (argc == 1)
    fname = argv[0];

  printf("[INFO] creating %s to hold the look-up tables for N = %ubits\n", fname, N);

  FILE *f = fopen(fname, "wb");
  if (f == NULL)
  {
    fprintf(stderr, "[ERROR] could not open file %s : %s", fname, strerror(errno));
    return 0;
  }

  submat R_i = {0};
  submat L_i = {0};
  submat I = {0};
  sqr_mat_I(&I);

  fwrite(&I, sizeof(I), 1, f);
  for (uint32_t i = 1; i < N; ++i)
  {
    sqr_mat_L_shift_n(&L_i, i);
    fwrite(&L_i, sizeof(L_i), 1, f);
  }

  fwrite(&I, sizeof(I), 1, f);
  for (uint32_t i = 1; i < N; ++i)
  {
    sqr_mat_R_shift_n(&R_i, i);
    fwrite(&R_i, sizeof(R_i), 1, f);
  }

  fclose(f);

  return 0;
}

#if N == 32U
#define _mm256_set1(uint32_t) _mm256_set1_epi32(uint32_t)
#define _mm256_setr(...) _mm256_setr_epi32(__VA_ARGS__)
#define _mm256_maskload(const_int_ptr, __m256i) _mm256_maskload_epi32((const int *)(const_int_ptr), __m256i)
#define _mm256_sub(__m256i1, __m256i2) _mm256_sub_epi32(__m256i1, __m256i2)
#define _mm256_sllv(__m256i1, __m256i2) _mm256_sllv_epi32(__m256i1, __m256i2)
#define _mm256_srli(__m256i, int) _mm256_srli_epi32(__m256i, int)
#define _mm256_slli(__m256i, int) _mm256_slli_epi32(__m256i, int)
#define SHIFTS_OFF_VALS 0, 1, 2, 3, 4, 5, 6, 7
#elif N == 64U
#define _mm256_set1(uint64_t) _mm256_set1_epi64x(uint64_t)
#define _mm256_setr(...) _mm256_setr_epi64x(__VA_ARGS__)
#define _mm256_maskload(const_uint64_ptr, __m256i) _mm256_maskload_epi64((const long long int *)(const_uint64_ptr), __m256i)
#define _mm256_sub(__m256i1, __m256i2) _mm256_sub_epi64(__m256i1, __m256i2)
#define _mm256_sllv(__m256i1, __m256i2) _mm256_sllv_epi64(__m256i1, __m256i2)
#define _mm256_srli(__m256i, int) _mm256_srli_epi64(__m256i, int)
#define _mm256_slli(__m256i, int) _mm256_slli_epi64(__m256i, int)
#define SHIFTS_OFF_VALS 0, 1, 2, 3
#else
#error "N must (currently) either 32 or 64"
#endif

// set ret to a square matrix that shifts right by n
int sqr_mat_R_shift_n(submat *ret, uint8_t n)
{
  __m256i shifts_off = _mm256_setr(SHIFTS_OFF_VALS);

  memset(ret, 0, sizeof(*ret));
  for (uint8_t i = 0; i * ROWS_PER_WIDTH < (N - n); ++i)
  {
    __m256i curr_shifts = _mm256_set1(N - n - (i * ROWS_PER_WIDTH) - 1);
    __m256i shifts = _mm256_sub(curr_shifts, shifts_off);
    MAT_ROWS_WIDTH(*ret, i) = _mm256_sllv(_mm256_set1(1), shifts);
  }

  return 0;
}

// set ret to a square matrix that shifts left by n
int sqr_mat_L_shift_n(submat *ret, uint8_t n)
{
  __m256i shifts_off = _mm256_setr(SHIFTS_OFF_VALS);

  memset(ret, 0, sizeof(*ret));
  for (uint8_t i = n / ROWS_PER_WIDTH; i * ROWS_PER_WIDTH < N; ++i)
  {
    __m256i curr_shifts = _mm256_set1(N - (i * ROWS_PER_WIDTH) - 1 + n);
    __m256i shifts = _mm256_sub(curr_shifts, shifts_off);
    MAT_ROWS_WIDTH(*ret, i) = _mm256_sllv(_mm256_set1(1), shifts);
  }

  return 0;
}

int sqr_mat_I(submat *ret)
{
  __m256i shifts_off = _mm256_setr(SHIFTS_OFF_VALS);

  memset(ret, 0, sizeof(*ret));
  for (uint8_t i = 0; i * ROWS_PER_WIDTH < N; ++i)
  {
    __m256i curr_shifts = _mm256_set1(N - (i * ROWS_PER_WIDTH) - 1);
    __m256i shifts = _mm256_sub(curr_shifts, shifts_off);
    MAT_ROWS_WIDTH(*ret, i) = _mm256_sllv(_mm256_set1(1), shifts);
  }

  return 0;
}

int sqr_mat_print(submat val)
{
  for (uint8_t j = 0; j < N; ++j)
  {
    for (uint8_t i = 0; i < N; ++i)
      printf("%c ", GET_MAT(val, i, j) ? '#' : '.');
    putchar('\n');
  }

  return 0;
}