#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <immintrin.h>
#include <windows.h>

#define __HELPER_IMPLEMENTATION__
#include "helper.h"

#define CANDIDATES_CAP (32 * 1024U)
#define SOLS_CAP (1024U)

#define FACTORS_CAP (256)

#define THREAD_COUNT (8U)

#define N (64U)

typedef __m256i mat_rows_t;
#if (N == 64U)
typedef uint64_t mat_row_t;
#define P (UINT64_MAX)
#elif (N == 32)
typedef uint32_t mat_row_t;
#define P ((uint64_t)UINT32_MAX)
#else
typedef uint64_t mat_row_t;
#define P (UINT64_MAX)
#error N must be (currently) either 32 or 64
#endif

#define BITS_PER_WIDTH ((uint32_t)sizeof(mat_rows_t) * CHAR_BIT)
#define ROWS_PER_WIDTH (BITS_PER_WIDTH / N)

typedef mat_rows_t sqr_mat[N / (ROWS_PER_WIDTH)];

#define MAT_ROWS_WIDTH(mat, i) ((mat)[(i)])

#define GET_BIT(dst, i) (((dst) >> (i)) & 0x1)
#define GET_MAT(mat, i, j) GET_BIT(((mat_row_t *)&(MAT_ROWS_WIDTH(mat, (mat_row_t)(j / ROWS_PER_WIDTH))))[j % ROWS_PER_WIDTH], N - (i)-1)

#define CLEAR_BIT(dst, i) ((dst) & (~(0x1 << (i))))

#define SET_BIT_TO_VAL(dst, i, val) \
  do                                \
  {                                 \
    (dst) = CLEAR_BIT((dst), (i));  \
    (dst) |= (val) << (i);          \
  } while (0)
// val MUST only have its lowest bit set
#define SET_MAT(mat, i, j, val) SET_BIT_TO_VAL(MAT_ROWS_WIDTH(mat, j), N - (i)-1, val)

typedef struct triplet_s
{
  uint8_t a, b, c;
} triplet;

typedef triplet candidates_list[CANDIDATES_CAP / (N - 1)];

typedef struct compute_sols_data_s
{
  candidates_list candidates;
  size_t used;
  uint16_t a;
} compute_sols_data_t;

int sqr_mat_L_shift_n(sqr_mat *ret, uint8_t n);
int sqr_mat_R_shift_n(sqr_mat *ret, uint8_t n);
int sqr_mat_I(sqr_mat *ret);
int sqr_mat_T(sqr_mat *res, uint32_t a, uint32_t b, uint32_t c);
int sqr_mat_print(sqr_mat val);
int sqr_mat_add(sqr_mat *res, sqr_mat a, sqr_mat b);
int sqr_mat_mult(sqr_mat *res, sqr_mat a, sqr_mat b);
int sqr_mat_pow_2_to_n(sqr_mat *res, sqr_mat a, uint32_t n); // computes a^(2^n)
int sqr_mat_bin_exp(sqr_mat *res, sqr_mat a, size_t n);      // computes a^n

int mat_row_print(mat_row_t row);
size_t prime_factors(uint32_t *factors, size_t arr_size, uint64_t n);

DWORD compute_sols_fixed_a(void *arg);
size_t compute_sols_fixed_a_inside(candidates_list candidates, uint16_t a);

int main(int argc, char **argv)
{
  (void)argc, (void)argv;

  timer global_time = {0};
  timer_start(&global_time);

  HANDLE thread_array[THREAD_COUNT] = {0};
  DWORD thread_id_array[THREAD_COUNT] = {0};

  // setup the data
  compute_sols_data_t candidates_datas[N - 1] = {0};
  for (uint16_t a = 1; a < N; ++a)
    candidates_datas[a - 1].a = a;

  for (uint16_t j = 0; j < (N - 1); j += THREAD_COUNT)
  {
    for (uint8_t i = 0; i < THREAD_COUNT; ++i)
    {
      if (i + j >= (N - 1))
        break;
      thread_array[i] = CreateThread(NULL, 0, &compute_sols_fixed_a, &(candidates_datas[j + i]), 0, &thread_id_array[i]);

      if (thread_array[i] == NULL)
      {
        fprintf(stderr, "[ERROR] while creating thread n %u!!", j + i);
        exit(69);
      }

      // tracks progress
      putchar('*');
    }
    WaitForMultipleObjects(THREAD_COUNT, thread_array, TRUE, INFINITE);
  }
  putchar('\n');

  // concatenate candidates into a single list
  size_t candidates_used_tot = 0;
  triplet candidates[CANDIDATES_CAP] = {0};
  for (uint16_t i = 0; i < N - 1; ++i)
  {
    memcpy(candidates + candidates_used_tot, candidates_datas[i].candidates, candidates_datas[i].used * sizeof(candidates[0]));
    candidates_used_tot += candidates_datas[i].used;
  }

  sqr_mat T = {0}, test = {0}, I = {0};
  sqr_mat_I(&I);

  // exclude candidates that have period < P
  triplet sols[SOLS_CAP] = {0};
  size_t sols_used = 0;

  uint32_t factors[FACTORS_CAP] = {0};
  const size_t num_prime_facts = prime_factors(factors, FACTORS_CAP, P);

  for (size_t i = 0; i < candidates_used_tot; ++i)
  {
    uint8_t valid = 1;
    for (size_t j = 0; j < num_prime_facts; ++j)
    {
      sqr_mat_T(&T, candidates[i].a, candidates[i].b, candidates[i].c);
      sqr_mat_bin_exp(&test, T, P / factors[j]);
      if (memcmp(&I, &test, sizeof(T)) == 0)
      {
        valid = 0;
        break;
      }
    }

    if (valid)
      sols[sols_used++] = candidates[i];
  }

  double tot_time = timer_stop(&global_time);
  printf("total time took %lf s\n", tot_time);

  printf("%llu sols found : ", sols_used);
  for (size_t i = 0; i < sols_used; ++i)
  {
    if (i % 9 == 0)
      printf("\n|");
    printf("%2u,%2u,%2u|", sols[i].a, sols[i].b, sols[i].c);
  }

  return 0;
}

DWORD compute_sols_fixed_a(void *arg)
{
  compute_sols_data_t *data = (compute_sols_data_t *)arg;
  data->used = compute_sols_fixed_a_inside(data->candidates, data->a);

  return 0;
}

size_t compute_sols_fixed_a_inside(candidates_list candidates, uint16_t a)
{
  size_t candidates_used = 0;

  sqr_mat T = {0};
  sqr_mat test = {0};

  sqr_mat I = {0};
  sqr_mat_I(&I);

  for (uint16_t b = 1; b < N; ++b)
  {
    for (uint16_t c = a; c < N; ++c)
    {
      sqr_mat_T(&T, a, b, c);

      sqr_mat_pow_2_to_n(&test, T, N);
      if (memcmp(&test, &T, sizeof(test)) == 0)
      {
        if (candidates_used > (CANDIDATES_CAP / (N - 1)))
        {
          // printf("%llu\n", candidates_used);
          fprintf(stderr, "[ERROR] buy mor RAM lol !!");
          return 69;
        }
        // assert(candidates_used < CANDIDATES_CAP / (N - 1));
        triplet sol = {.a = a, .b = b, .c = c};
        candidates[candidates_used++] = sol;
      }
    }
  }

  // printf("%llu candidates found : ", candidates_used);
  // for (size_t i = 0; i < candidates_used; ++i)
  //   printf("%u, %u, %u | ", candidates[i].a, candidates[i].b, candidates[i].c);
  // putchar('\n');

  return candidates_used;
}

#if N == 64
#define _mm256_set1(uint64_t) _mm256_set1_epi64x(uint64_t)
#define _mm256_setr(...) _mm256_setr_epi64x(__VA_ARGS__)
#define _mm256_maskload(const_uint64_ptr, __m256i) _mm256_maskload_epi64((const long long int *)(const_uint64_ptr), __m256i)
#define _mm256_sub(__m256i1, __m256i2) _mm256_sub_epi64(__m256i1, __m256i2)
#define _mm256_sllv(__m256i1, __m256i2) _mm256_sllv_epi64(__m256i1, __m256i2)
#define _mm256_srli(__m256i, int) _mm256_srli_epi64(__m256i, int)
#define _mm256_slli(__m256i, int) _mm256_slli_epi64(__m256i, int)
#define SHIFTS_OFF_VALS 0, 1, 2, 3
#else
#define _mm256_set1(uint32_t) _mm256_set1_epi32(uint32_t)
#define _mm256_setr(...) _mm256_setr_epi32(__VA_ARGS__)
#define _mm256_maskload(const_int_ptr, __m256i) _mm256_maskload_epi32((const int *)(const_int_ptr), __m256i)
#define _mm256_sub(__m256i1, __m256i2) _mm256_sub_epi32(__m256i1, __m256i2)
#define _mm256_sllv(__m256i1, __m256i2) _mm256_sllv_epi32(__m256i1, __m256i2)
#define _mm256_srli(__m256i, int) _mm256_srli_epi32(__m256i, int)
#define _mm256_slli(__m256i, int) _mm256_slli_epi32(__m256i, int)
#define SHIFTS_OFF_VALS 0, 1, 2, 3, 4, 5, 6, 7
#endif // N == 64

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

int sqr_mat_T(sqr_mat *res, uint32_t a, uint32_t b, uint32_t c)
{
  sqr_mat L_a = {0};
  sqr_mat R_b = {0};
  sqr_mat L_c = {0};
  sqr_mat I = {0};
  sqr_mat T1 = {0}, T2 = {0}, T3 = {0};
  sqr_mat T_ = {0};

  sqr_mat_I(&I);

  sqr_mat_L_shift_n(&L_a, a);
  sqr_mat_L_shift_n(&L_c, c);
  sqr_mat_R_shift_n(&R_b, b);
  sqr_mat_I(&I);

  sqr_mat_add(&T1, I, L_a);
  sqr_mat_add(&T2, I, R_b);
  sqr_mat_add(&T3, I, L_c);

  sqr_mat_mult(&T_, T1, T2);
  sqr_mat_mult(res, T_, T3);

  return 0;
}

int sqr_mat_add(sqr_mat *res, sqr_mat a, sqr_mat b)
{
  for (uint8_t i = 0; i < (N / ROWS_PER_WIDTH); ++i)
    MAT_ROWS_WIDTH(*res, i) = _mm256_xor_si256(MAT_ROWS_WIDTH(a, i), MAT_ROWS_WIDTH(b, i));
  return 0;
}

int sqr_mat_mult(sqr_mat *res, sqr_mat a, sqr_mat b)
{
  __m256i mask = _mm256_set1(((mat_row_t)1) << (N - 1)); // 0b100000...000

  for (uint32_t i = 0; i < (N / ROWS_PER_WIDTH); ++i)
  {
    MAT_ROWS_WIDTH(*res, i) = _mm256_set1(0);
    for (uint32_t j = 0; j < N; ++j)
    {
      __m256i mask_a = _mm256_and_si256(_mm256_slli(MAT_ROWS_WIDTH(a, i), j), mask); // the highest bit tells if the row from b should be xored with res

      __m256i v_b_j = _mm256_set1(((mat_row_t *)b)[j]);

      __m256i vals = _mm256_maskload((mat_row_t *)&v_b_j, mask_a);
      MAT_ROWS_WIDTH(*res, i) = _mm256_xor_si256(MAT_ROWS_WIDTH(*res, i), vals);
    }
  }
  return 0;
}

int sqr_mat_pow_2_to_n(sqr_mat *res, sqr_mat a, uint32_t n)
{
  sqr_mat tmp = {0};
  memcpy(&tmp, a, sizeof(tmp));
  for (uint32_t i = 0; i < n; ++i)
  {
    sqr_mat_mult(res, tmp, tmp);
    memcpy(&tmp, res, sizeof(tmp));
  }
  return 0;
}

int sqr_mat_bin_exp(sqr_mat *res, sqr_mat a, size_t n)
{
  sqr_mat_I(res);

  sqr_mat a_ = {0};
  memcpy(&a_, a, sizeof(a_));
  while (n > 0)
  {
    if (n & 1)
    {
      // res = res * a
      sqr_mat tmp = {0};
      sqr_mat_mult(&tmp, *res, a_);
      memcpy(res, &tmp, sizeof(*res));
    }
    // a = a * a;
    sqr_mat tmp = {0};
    sqr_mat_mult(&tmp, a_, a_);
    memcpy(&a_, &tmp, sizeof(a_));

    n >>= 1;
  }

  return 0;
}

int sqr_mat_print(sqr_mat val)
{
  for (uint8_t j = 0; j < N; ++j)
  {
    for (uint8_t i = 0; i < N; ++i)
      printf("%c ", GET_MAT(val, i, j) ? '#' : '.');
    putchar('\n');
  }

  return 0;
}

#undef GET_BIT
#define GET_BIT(val, idx) ((val) >> (idx)) & 0x1

int mat_row_print(mat_row_t row)
{
  for (uint8_t i = 0; i < (sizeof(row) / sizeof(uint8_t)); ++i)
  {
    uint8_t bits = row >> ((sizeof(row) * CHAR_BIT) - ((i + 1) * 8));
    for (uint8_t h = 0; h < 2; ++h)
    {
      for (uint8_t j = h * 4; j < (h + 1) * 4; ++j)
        printf("%u", GET_BIT(bits, CHAR_BIT - (j + 1)));
      putchar(' ');
    }
  }
  return 0;
}

size_t prime_factors(uint32_t *factors, size_t arr_size, uint64_t n)
{
  assert(arr_size > 0);
  size_t arr_used = 0;

  if (n % 2 == 0)
  {
    assert(arr_size > arr_used);
    factors[arr_used++] = 2;
    n /= 2;
    while (n % 2 == 0)
      n /= 2;
  }

  for (int i = 3; i * i <= n; i = i + 2)
  {
    if (n % i == 0)
    {
      assert(arr_size > arr_used);
      factors[arr_used++] = i;
      n /= i;
      while (n % i == 0)
        n /= i;
    }
  }

  if (n > 2)
  {
    assert(arr_size > arr_used);
    factors[arr_used++] = n;
  }

  return arr_used;
}