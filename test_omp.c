// gcc -fopenmp -O2 test_omp.c -Lutility -lutility -o test_omp
// export OMP_NUM_THREADS=4 ; ./test_omp

#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "headers/sort_algs.h"
#include "headers/test_gen.h"
#include "headers/test_sizes.h"

#define ITERATIONS_PER_SIZE 10
#define MIN_SIZE 10
#define MID_SIZE 10
// #define MAX_SIZE 100000000
#define MAX_SIZE 100
#define MULTIPLIER 1.5

#define SWAP(type, temp, x, y) type(temp) = (x); (x) = (y); (y) = (temp)

ValueType *MergeSortOMP(ValueType *array, ValueType *auxil, const IndexType n) {
  for (IndexType sort_len = 2; sort_len < (n << 1); sort_len <<= 1) {
#pragma omp parallel for
    for (IndexType start1 = 0; start1 < n; start1 += sort_len) {
      IndexType start2 = start1 + (sort_len >> 1);
      IndexType end1 = start2, end2 = start1 + sort_len, i1 = start1, i2 = start2, j = start1;
      if (end1 > n) {
        end1 = n;
      }
      if (end2 > n) {
        end2 = n;
      }
      if (i2 < end2) {
        while (1) {
          if (array[i2] < array[i1]) {
            auxil[j++] = array[i2++];
            if (i2 == end2) {
              break;
            }
          } else {
            auxil[j++] = array[i1++];
            if (i1 == end1) {
              break;
            }
          }
        }
      }
      while (i1 < end1) {
        auxil[j++] = array[i1++];
      }
      while (i2 < end2) {
        auxil[j++] = array[i2++];
      }
    }
    SWAP(ValueType *, temp, array, auxil);
  }
  return array;
}

int main() {
  IndexType num_sizes;
  IndexType* sizes = GetTestSizes(MIN_SIZE, MID_SIZE, MAX_SIZE, MULTIPLIER, &num_sizes);

  // TEST SIZE LOOP
  for (IndexType itsize = 0; itsize < num_sizes; ++itsize) {
    IndexType n = sizes[itsize];
    ValueType* original = GetSortedInts(n);
    ValueType* test = (ValueType*)malloc(n * sizeof(ValueType));
    ValueType* buff = (ValueType*)malloc(n * sizeof(ValueType));
    double average = 0.0, variance = 0.0;
    double times[ITERATIONS_PER_SIZE];

    // ITERATION LOOP
    for (IndexType iter = 0; iter < ITERATIONS_PER_SIZE; ++iter) {
      memcpy(test, original, n * sizeof(ValueType));
      test = GetRandomizeTest(test, n, 3 + itsize + iter);
      double t = omp_get_wtime();

      // MAIN CALCULATION
      ValueType* res = MergeSortOMP(test, buff, n);

      t = omp_get_wtime() - t;
      average += t;
      times[iter] = t;
    }

    // TIME CALCULATION
    average /= ITERATIONS_PER_SIZE;
    if (ITERATIONS_PER_SIZE > 1) {
      for (IndexType i = 0; i < ITERATIONS_PER_SIZE; ++i) {
        times[i] -= average;
        variance += times[i] * times[i];
      }
      variance = MySqrt(variance / (ITERATIONS_PER_SIZE - 1));
    }

    // PRINT RESULT FOR SIZE
    printf("OMP MergeSort: N=%d num_proc=%d max_th=%d time=%lf +- %lf\n", n, omp_get_num_procs(),
         omp_get_max_threads(), average, variance);

    free(buff);
    free(test);
    free(original);
  }

  free(sizes);
  return 0;
}