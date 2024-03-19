#define __USE_MINGW_ANSI_STDIO 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <immintrin.h>
#include <windows.h>

#define __TIMER_IMPLEMENTATION__
#include "timer.h"

#include "constants.h"

#define __SQUARE_MATRIX_ENABLE_WARNINGS__

#if N > 64
#define __SQUARE_MATRIX_IMPLEMENTATION__
#define __SQUARE_MATRIX_SIZE__ (SUBN)
#define __SQUARE_MATRIX_TYPE__ submat
#define __SQUARE_MATRIX_ROW_TYPE__ subrow_t
#define __SQUARE_MATRIX_ROWS_PER_WIDTH__ SUBROWS_PER_SUBWIDTH
#define __SQUARE_MATRIX_GET_MAT__ GET_SUBMAT
#include "square_matrix.h"
#endif

#define __SQUARE_MATRIX_IMPLEMENTATION__
#define __SQUARE_MATRIX_SIZE__ (N)
#define __SQUARE_MATRIX_SUBROW_COUNT__ SUBROWS_PER_ROW
#define __SQUARE_MATRIX_TYPE__ sqr_mat
#define __SQUARE_MATRIX_ROW_TYPE__ mat_row_t
#define __SQUARE_MATRIX_SUBMAT_NAMESPACE__ submat_
#define __SQUARE_MATRIX_ROWS_PER_WIDTH__ ROWS_PER_WIDTH
#define __SQUARE_MATRIX_GET_MAT__ GET_MAT
#ifdef USE_SUBMAT_MULT
#define __SQUARE_MATRIX_MULT_SUBMAT__
#endif
#include "square_matrix.h"

#define CANDIDATES_CAP (32 * 1024U)
#define SOLS_CAP (1024U)

#define FACTORS_CAP (256)

#define THREAD_COUNT (SUBN - 1)

#define LOG_SUCC_P (N)

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

#if N == 32U
size_t num_facts_of_P = 5;
uint64_t factors_of_P[] = {3, 5, 17, 257, 65537};

#elif N == 64U
size_t num_facts_of_P = 7;
uint64_t factors_of_P[] = {3, 5, 17, 257, 641, 65537, 6700417};

#elif N == 128U
size_t num_facts_of_P = 9;
uint64_t factors_of_P[] = {3, 5, 17, 257, 641, 65537, 274177, 6700417, 67280421310721};
#else
#error "N must be (currently) either 32, 64 or 128"
#endif

submat R_i_lut[SUBN];
submat L_i_lut[SUBN];

int sqr_mat_T(sqr_mat *res, uint32_t a, uint32_t b, uint32_t c);
int is_2_to_N_periodic(sqr_mat *T);
int submat_copy(sqr_mat *dst, submat src, size_t pos);
int load_lookup_tables(submat (*L_i)[SUBN], submat (*R_i)[SUBN], char *fname);

DWORD compute_sols_fixed_a(void *arg);
size_t compute_sols_fixed_a_inside(sols_list candidates, uint16_t a);

int main(int argc, char **argv)
{
  --argc, ++argv; // remove prog name

  /* -------------------- loading the look-up table ----------------------*/

  // default name for the look-up table "lut_xx.data"
  char lut_fname[256] = LUT_FNAME;
  if (argc > 0)
    strncpy(lut_fname, argv[0], 255);

  if (load_lookup_tables(&L_i_lut, &R_i_lut, lut_fname))
    return 69;

  // comfirmation line
  printf("[INFO] computing xorshift coefficients for N = %ubits implemented using %" PRIu64 " uint%u_t%c, with period, P = 2^%u - 1, using %u threads\n", N, SUBROWS_PER_ROW, SUBN, SUBROWS_PER_ROW > 1 ? 's' : 0, N, THREAD_COUNT);

  // setup and start timer
  timer global_time = {0};
  timer_start(&global_time);

  /* ------------ prepare for multithreading ----------------*/

  HANDLE thread_array[THREAD_COUNT] = {0};
  DWORD thread_id_array[THREAD_COUNT] = {0};

  // setup the data
  compute_sols_data_t sols_datas[SUBN - 1] = {0};
  for (uint16_t a = 1; a < SUBN; ++a)
    sols_datas[a - 1].a = a;

  // start thread by batches of THREAD_COUNT
  for (uint16_t j = 0; j < (SUBN - 1); j += THREAD_COUNT)
  {
    for (uint8_t i = 0; i < THREAD_COUNT; ++i)
    {
      if (i + j >= (SUBN - 1))
        break;

      // each thread tries `b` and `c` for its own `a`
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
  for (uint16_t i = 0; i < SUBN - 1; ++i)
  {
    memcpy(sols + sols_used_tot, sols_datas[i].sols, sols_datas[i].used * sizeof(sols[0]));
    sols_used_tot += sols_datas[i].used;
  }

  // stop timer and report time
  double tot_time = timer_stop(&global_time);
  printf("total time took %lf s\n", tot_time);

  /* ----------------- reporting the data ----------------------*/

  printf("found %" PRIu64 " triplets such that T = (I + L^a)(I + R^b)(I + L^b) has order of 2^%u - 1 = %" PRIu64 ": ", sols_used_tot, N, P);
  // printing the solutions by groups of nine
  for (size_t i = 0; i < sols_used_tot; ++i)
  {
    if (i % 9 == 0)
      printf("\n|");
    printf("%2u,%2u,%2u|", sols[i].a, sols[i].b, sols[i].c);
  }

  return 0;
}

// wrapper for win32 API
DWORD compute_sols_fixed_a(void *arg)
{
  compute_sols_data_t *data = (compute_sols_data_t *)arg;
  data->used = compute_sols_fixed_a_inside(data->sols, data->a);

  return 0;
}

/* finds values of b and c such that T(a, b, c) has period 2^N - 1 */
size_t compute_sols_fixed_a_inside(sols_list sols, uint16_t a)
{
  size_t sols_used = 0;

  sqr_mat T = {0};

  for (uint16_t b = 1; b < SUBN; ++b)
  {
#if N > 64
    for (uint16_t c = 1; c < SUBN; ++c)
#else
    for (uint16_t c = a + 1; c < SUBN; ++c)
#endif
    {
      sqr_mat_T(&T, a, b, c); // T = (I + L^a)(I + R^b)(I + L^c)

      if (is_2_to_N_periodic(&T))
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

  return sols_used;
}

int is_2_to_N_periodic(sqr_mat *T)
{
  sqr_mat test = {0};
#if N > 64
  sqr_mat tmp = {0};
#endif

  sqr_mat_pow_2_to_n(&test, *T, LOG_SUCC_P);
  if (memcmp(&test, T, sizeof(sqr_mat)) != 0)
    return 0;

  // exclude candidates that have period < P
  for (size_t j = 0; j < num_facts_of_P; ++j)
  {
#if N > 64
    memcpy(&tmp, T, sizeof(sqr_mat));
    for (uint32_t i = 0; i < num_facts_of_P; ++i)
    {
      if (i == j)
        continue;
      sqr_mat_bin_exp(&test, tmp, factors_of_P[i]);
      memcpy(T, &test, sizeof(sqr_mat));
    }
#else
    sqr_mat_bin_exp(&test, *T, P / factors_of_P[j]);
#endif
    if (memcmp(&I, &test, sizeof(sqr_mat)) == 0)
      return 0;
  }

  return 1;
}

int sqr_mat_T(sqr_mat *res, uint32_t a, uint32_t b, uint32_t c)
{
  submat T1 = {0};
  submat T2 = {0};
  submat T3 = {0};
  submat T_ = {0};

#if N <= 64U
  sqr_mat_add(&T1, I, L_i_lut[a]);
  sqr_mat_add(&T2, I, R_i_lut[b]);
  sqr_mat_add(&T3, I, L_i_lut[c]);

  sqr_mat_mult(&T_, T1, T2);
  sqr_mat_mult(res, T_, T3);
#else
  submat_add(&T1, I, L_i_lut[a]);
  submat_add(&T2, I, R_i_lut[b]);
  submat_add(&T3, I, R_i_lut[c]);

  submat_mult(&T_, T1, T2);

  submat_copy(res, T_, SUBROWS_PER_ROW - 1);
  submat_copy(res, T3, (SUBROWS_PER_ROW * SUBROWS_PER_ROW - 1));

  for (uint32_t i = 1; i < SUBROWS_PER_ROW; ++i)
    submat_copy(res, I, i * SUBROWS_PER_ROW + (i - 1));
#endif

  return 0;
}

#if N > 64
int submat_copy(sqr_mat *dst, submat src, size_t pos)
{
  assert(pos < (SUBROWS_PER_ROW * SUBROWS_PER_ROW));

#ifndef USE_SUBMAT_MULT
  subrow_t src_width[SUBROWS_PER_SUBWIDTH];
  // submat_print(src);
  // printf("\n---------------------------------\n\n");
  const size_t off = pos % SUBROWS_PER_ROW + ((pos / SUBROWS_PER_ROW) * SUBROWS_PER_ROW * SUBN);
  subrow_t *mat = (subrow_t *)*dst;
  for (uint32_t j = 0; j < SUBN / SUBROWS_PER_SUBWIDTH; ++j)
  {
    _mm256_storeu_si256((__m256i_u *)src_width, MAT_ROWS_WIDTH(src, j));
    for (uint32_t i = 0; i < SUBROWS_PER_SUBWIDTH; ++i)
      mat[off + (j * SUBROWS_PER_WIDTH + i) * SUBROWS_PER_ROW] = src_width[i];
  }
#else
  memcpy(&((*dst)[pos]), src, sizeof((*dst)[pos]));
#endif

  return 0;
}
#endif

/* parentheses are necessary: https://cdecl.org/?q=int+%28*R_i%29%5B32%5D */
int load_lookup_tables(submat (*L_i)[SUBN], submat (*R_i)[SUBN], char *fname)
{
  FILE *f_L = fopen(fname, "rb");
  FILE *f_R = fopen(fname, "rb");
  if (f_L == NULL || f_R == NULL)
  {
    fprintf(stderr, "something went wrong when opening \"%s\" : %s\n", fname, strerror(errno));
    return 1;
  }
  fseek(f_R, SUBN * sizeof(submat), SEEK_SET); // jump to R_i data

  // loading the look-up tables
  fread(L_i, sizeof(submat), SUBN, f_L);
  fread(R_i, sizeof(submat), SUBN, f_R);

  fclose(f_L);
  fclose(f_R);

  return 0;
}