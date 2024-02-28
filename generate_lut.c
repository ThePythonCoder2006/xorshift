#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "constants.h"

#undef I

int sqr_mat_R_shift_n(sqr_mat *ret, uint8_t n);
int sqr_mat_L_shift_n(sqr_mat *ret, uint8_t n);
int sqr_mat_I(sqr_mat *ret);
int write_sqr_mat(FILE *f, sqr_mat mat);

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    fprintf(stderr, "[USAGE] %s <output_file_txt>\n", argv[0]);
    return 1;
  }

  --argc, ++argv; // remove the program name

  FILE *f = fopen(argv[0], "w");
  if (f == NULL)
  {
    fprintf(stderr, "[ERROR] could not open file %s : %s", argv[0], strerror(errno));
    return 0;
  }

  sqr_mat R_i = {0};
  sqr_mat L_i = {0};
  sqr_mat I = {0};
  sqr_mat_I(&I);

  // header
  fprintf(f, "#include \"constants.h\"\n\nsqr_mat R_i_lut[N] = {\n");
  write_sqr_mat(f, I);
  fprintf(f, ",\n");
  for (uint32_t i = 1; i < N; ++i)
  {
    sqr_mat_R_shift_n(&R_i, i);
    write_sqr_mat(f, R_i);
    fprintf(f, ",\n");
  }
  fprintf(f, "};");

  fprintf(f, "\n\nsqr_mat L_i_lut[N] = {\n");
  write_sqr_mat(f, I);
  fprintf(f, ",\n");
  for (uint32_t i = 1; i < N; ++i)
  {
    sqr_mat_L_shift_n(&L_i, i);
    write_sqr_mat(f, L_i);
    fprintf(f, ",\n");
  }
  fprintf(f, "};");

  fclose(f);

  return 0;
}

#if N == 32
#define _mm256_set1(uint32_t) _mm256_set1_epi32(uint32_t)
#define _mm256_setr(...) _mm256_setr_epi32(__VA_ARGS__)
#define _mm256_maskload(const_int_ptr, __m256i) _mm256_maskload_epi32((const int *)(const_int_ptr), __m256i)
#define _mm256_sub(__m256i1, __m256i2) _mm256_sub_epi32(__m256i1, __m256i2)
#define _mm256_sllv(__m256i1, __m256i2) _mm256_sllv_epi32(__m256i1, __m256i2)
#define _mm256_srli(__m256i, int) _mm256_srli_epi32(__m256i, int)
#define _mm256_slli(__m256i, int) _mm256_slli_epi32(__m256i, int)
#define SHIFTS_OFF_VALS 0, 1, 2, 3, 4, 5, 6, 7
#elif N == 64
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
int sqr_mat_R_shift_n(sqr_mat *ret, uint8_t n)
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
int sqr_mat_L_shift_n(sqr_mat *ret, uint8_t n)
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

int sqr_mat_I(sqr_mat *ret)
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

// writes the contents of a square matrix into a c file in a format making it includable by the main program
int write_sqr_mat(FILE *f, sqr_mat mat)
{
  putc('{', f);
  for (uint8_t i = 0; i < N / ROWS_PER_WIDTH; ++i)
  {
    putc('{', f);
    for (uint8_t j = 0; j < sizeof(mat_width_t) / sizeof(uint64_t); ++j)
    {
      fprintf(f, "0x%016llxULL", ((uint64_t *)&(MAT_ROWS_WIDTH(mat, i)))[j]);

      if (j != (sizeof(mat_width_t) / sizeof(uint64_t)) - 1)
        fprintf(f, ", ");
    }
    putc('}', f);

    // if (i != (N / ROWS_PER_WIDTH) - 1)
    fprintf(f, ",\n");
  }
  putc('}', f);

  return 0;
}