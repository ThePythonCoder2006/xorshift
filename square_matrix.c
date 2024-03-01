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

#include "constants.h"

#define CANDIDATES_CAP (32 * 1024U)
#define SOLS_CAP (1024U)

#define FACTORS_CAP (256)

#define THREAD_COUNT (N - 1)

#define LOG_SUCC_P (N)

#define GET_BIT(dst, i) (((dst) >> (i)) & 0x1)
#if N == 128U
#define GET_MAT(mat, i, j) GET_BIT(((uint64_t *)&(MAT_ROWS_WIDTH(mat, (uint64_t)((j) / (2 * ROWS_PER_WIDTH)) + ((j) & 1))))[j % ROWS_PER_WIDTH], N - (i)-1)
#else
#define GET_MAT(mat, i, j) GET_BIT(((mat_row_t *)&(MAT_ROWS_WIDTH(mat, (uint64_t)(j / (ROWS_PER_WIDTH)))))[j % ROWS_PER_WIDTH], N - (i)-1)
#endif

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

typedef triplet sols_list[CANDIDATES_CAP / (N - 1)];

typedef struct compute_sols_data_s
{
  sols_list sols;
  size_t used;
  uint16_t a;
} compute_sols_data_t;

int sqr_mat_T(sqr_mat *res, uint32_t a, uint32_t b, uint32_t c);
int sqr_mat_print(sqr_mat val);
int sqr_mat_add(sqr_mat *res, sqr_mat a, sqr_mat b);
int sqr_mat_mult(sqr_mat *res, sqr_mat a, sqr_mat b);
int sqr_mat_pow_2_to_n(sqr_mat *res, sqr_mat a, uint32_t n); // computes a^(2^n)
int sqr_mat_bin_exp(sqr_mat *res, sqr_mat a, size_t n);      // computes a^n

// int mat_row_print(mat_row_t row);
size_t prime_factors(uint32_t *factors, size_t arr_size, uint64_t n);

DWORD compute_sols_fixed_a(void *arg);
size_t compute_sols_fixed_a_inside(sols_list candidates, uint16_t a);

int main(int argc, char **argv)
{
  (void)argc, (void)argv;

  sqr_mat_print(I);

  return 0;

  // comfirmation line
  printf("[INFO] computing xorshift coefficients for N = %u => P = 0x%016llx and LOG_SUCC_P = %u, using %u threads\n", N, P, LOG_SUCC_P, THREAD_COUNT);

  timer global_time = {0};
  timer_start(&global_time);

  HANDLE thread_array[THREAD_COUNT] = {0};
  DWORD thread_id_array[THREAD_COUNT] = {0};

  // setup the data
  compute_sols_data_t sols_datas[N - 1] = {0};
  for (uint16_t a = 1; a < N; ++a)
    sols_datas[a - 1].a = a;

  // main loops
  for (uint16_t j = 0; j < (N - 1); j += THREAD_COUNT)
  {
    for (uint8_t i = 0; i < THREAD_COUNT; ++i)
    {
      if (i + j >= (N - 1))
        break;

      // each thread tries b and c for its own a
      thread_array[i] = CreateThread(NULL, 0, &compute_sols_fixed_a, &(sols_datas[j + i]), 0, &thread_id_array[i]);

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
  size_t sols_used_tot = 0;
  triplet sols[CANDIDATES_CAP] = {0};
  for (uint16_t i = 0; i < N - 1; ++i)
  {
    memcpy(sols + sols_used_tot, sols_datas[i].sols, sols_datas[i].used * sizeof(sols[0]));
    sols_used_tot += sols_datas[i].used;
  }

  double tot_time = timer_stop(&global_time);
  printf("total time took %lf s\n", tot_time);

  printf("found %llu triplets such that T = (I + L^a)(I + R^b)(I + L^b) has order of 2^%u - 1 = %llu: ", sols_used_tot, N, P);
  // printing the solutions by groups of nine
  for (size_t i = 0; i < sols_used_tot; ++i)
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
  data->used = compute_sols_fixed_a_inside(data->sols, data->a);

  return 0;
}

size_t compute_sols_fixed_a_inside(sols_list sols, uint16_t a)
{
  size_t sols_used = 0;

  sqr_mat T = {0};
  sqr_mat test = {0};

  for (uint16_t b = 1; b < N; ++b)
  {
    for (uint16_t c = a; c < N; ++c)
    {
      sqr_mat_T(&T, a, b, c);

      sqr_mat_pow_2_to_n(&test, T, LOG_SUCC_P);
      if (memcmp(&test, &T, sizeof(test)) == 0)
      {
        // exclude candidates that have period < P
        uint8_t valid = 1;
        for (size_t j = 0; j < num_facts_of_P; ++j)
        {
          sqr_mat_T(&T, a, b, c); // T = (I + L^a)(I + R^b)(I + L^c)
          sqr_mat_bin_exp(&test, T, P / factors_of_P[j]);
          if (memcmp(&I, &test, sizeof(T)) == 0)
          {
            valid = 0;
            break;
          }
        }

        if (valid)
        {
          if (sols_used > (CANDIDATES_CAP / (N - 1)))
          {
            fprintf(stderr, "[ERROR] buy mor RAM lol !!");
            return 69;
          }

          sols[sols_used++] = (triplet){.a = a, .b = b, .c = c};
        }
      }
    }
  }

  // printf("%llu candidates found : ", candidates_used);
  // for (size_t i = 0; i < candidates_used; ++i)
  //   printf("%u, %u, %u | ", candidates[i].a, candidates[i].b, candidates[i].c);
  // putchar('\n');

  return sols_used;
}

#define _mm256_set_zero() _mm256_set1_epi64x(0)
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
#define _mm256_maskload(const_int_ptr, __m256i) _mm256_maskload_epi64((const long long int *)(const_int_ptr), __m256i)
#define _mm256_sub(__m256i1, __m256i2) _mm256_sub_epi64(__m256i1, __m256i2)
#define _mm256_sllv(__m256i1, __m256i2) _mm256_sllv_epi64(__m256i1, __m256i2)
#define _mm256_srli(__m256i, int) _mm256_srli_epi64(__m256i, int)
#define _mm256_slli(__m256i, int) _mm256_slli_epi64(__m256i, int)
#define SHIFTS_OFF_VALS 0, 1, 2, 3
#elif N == 128
#define _mm256_slli(__m256i, int) _mm256_slli_epi16(_mm256_bslli_epi128(__m256i, int >> 3), int & 7) /* int >> 3 = int / 8 and int & 7 == int % 8*/
#define _mm256_maskload(const_int_ptr, __m256i_mask) _mm256_maskload_epi64((const long long int *)(const_int_ptr), _mm256_permute4x64_epi64(__m256i_mask, (uint8_t)0xA5))
#define _mm256_set1(__m128i) _mm256_broadcastsi128_si256(__m128i)
#else
#error "N must (currently) either 32 or 64"
#endif

int sqr_mat_T(sqr_mat *res, uint32_t a, uint32_t b, uint32_t c)
{
  sqr_mat L_a;
  memcpy(&L_a, &L_i_lut[a], sizeof(sqr_mat));
  sqr_mat R_b;
  memcpy(&R_b, &R_i_lut[b], sizeof(sqr_mat));
  sqr_mat L_c;
  memcpy(&L_c, &L_i_lut[c], sizeof(sqr_mat));
  sqr_mat T1 = {0};
  sqr_mat T2 = {0};
  sqr_mat T3 = {0};
  sqr_mat T_ = {0};

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
#if N == 128U
  __m256i mask = _mm256_set_epi64x(0x8000000000000000ULL, 0ULL, 0x8000000000000000ULL, 0ULL);
#else
  __m256i mask = _mm256_set1(((mat_row_t)1) << (N - 1)); // 0b100000...000
#endif

  for (uint32_t i = 0; i < (N / ROWS_PER_WIDTH); ++i)
  {
    MAT_ROWS_WIDTH(*res, i) = _mm256_set_zero();

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
  memcpy(res, &I, sizeof(*res));

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

// #undef GET_BIT
// #define GET_BIT(val, idx) ((val) >> (idx)) & 0x1

// int mat_row_print(mat_row_t row)
// {
//   for (uint8_t i = 0; i < (sizeof(row) / sizeof(uint8_t)); ++i)
//   {
//     uint8_t bits = row >> ((sizeof(row) * CHAR_BIT) - ((i + 1) * 8));
//     for (uint8_t h = 0; h < 2; ++h)
//     {
//       for (uint8_t j = h * 4; j < (h + 1) * 4; ++j)
//         printf("%u", GET_BIT(bits, CHAR_BIT - (j + 1)));
//       putchar(' ');
//     }
//   }
//   return 0;
// }

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