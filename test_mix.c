// mpicc -fopenmp -O2 test_mix.c -Lutility -lutility -o test_mix
// export OMP_NUM_THREADS=4 ; mpirun -np 4 ./test_mix

#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "headers/ordering.h"
#include "headers/sort_algs.h"
#include "headers/test_gen.h"
#include "headers/test_sizes.h"

#define ITERATIONS_PER_SIZE 10
#define MIN_SIZE 10
#define MID_SIZE 10
#define MAX_SIZE 100000000
// #define MAX_SIZE 100
#define MULTIPLIER 1.3

#define SWAP(type, temp, x, y) \
  type(temp) = (x);            \
  (x) = (y);                   \
  (y) = (temp)

ValueType* MergeSortOMP(ValueType* array, ValueType* auxil, const IndexType n) {
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
    SWAP(ValueType*, temp, array, auxil);
  }
  return array;
}

int main(int argc, char** argv) {
  int id, np;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  IndexType num_sizes;
  IndexType* sizes = GetTestSizes(MIN_SIZE, MID_SIZE, MAX_SIZE, MULTIPLIER, &num_sizes);

  SortNet ordering;
  struct MyVector order;

  ///////////////////////////////////////////
  ////  O R D E R I N G     C H O I C E  ////
  ////      0   :   O D D - E V E N      ////
  ////      1   :   B I T O N I C        ////
  ///////////////////////////////////////////

  for (short current_order = 0; current_order < 2; ++current_order) {
    ordering = current_order ? BitonicOrdering(np) : OddEvenOrdering(np);
    order = ordering[id];
    for (short i = 0; i < np; ++i) {
      if (i != id) {
        free(ordering[i].data);
      }
    }
    free(ordering);

    // TEST SIZE LOOP
    for (IndexType itsize = 0; itsize < num_sizes; ++itsize) {
      IndexType n_all = sizes[itsize];
      IndexType n_one = (IndexType)((n_all - 1) / np) + 1;
      n_all = np * n_one;

      ValueType* test = (ValueType*)malloc(n_one * sizeof(ValueType));
      ValueType* buff = (ValueType*)malloc(n_one * sizeof(ValueType));
      ValueType* recv = (ValueType*)malloc(n_one * sizeof(ValueType));

      double average = 0.0, variance = 0.0;
      double times[ITERATIONS_PER_SIZE];

      // ITERATION LOOP
      for (IndexType iter = 0; iter < ITERATIONS_PER_SIZE; ++iter) {
        SetRandomInts(test, n_one, -n_all, n_all, 3 + n_one + iter + ITERATIONS_PER_SIZE * id);
        MPI_Barrier(MPI_COMM_WORLD);
        double t = MPI_Wtime();

        // MAIN CALCULATION
        // Sort Own
        ValueType* res = MergeSortOMP(test, buff, n_one);
        if (res != test) {
          SWAP(ValueType*, temp, test, buff);
        }
        // Merge Sorted
        for (short iord = 0; iord < order.size; ++iord) {
          MPI_Sendrecv(test, n_one, MPI_INT, order.data[iord], 0, recv, n_one, MPI_INT, order.data[iord], 0,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          if (id < order.data[iord]) {
            IndexType i = 0, j = 0;
            for (IndexType k = 0; k < n_one; ++k) {
              if (recv[j] < test[i]) {
                buff[k] = recv[j++];
              } else {
                buff[k] = test[i++];
              }
            }
          } else {
            IndexType i = n_one - 1, j = n_one - 1;
            for (IndexType k = n_one - 1; k != -1; --k) {
              if (test[i] < recv[j]) {
                buff[k] = recv[j--];
              } else {
                buff[k] = test[i--];
              }
            }
          }
          SWAP(ValueType*, temp, test, buff);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        t = MPI_Wtime() - t;

        if (id == 0) {
          average += t;
          times[iter] = t;
        }
      }

      if (id == 0) {
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
        if (current_order) {
          printf("MIX Bitonic MergeSort: N=%d np=%d num_proc=%d max_th=%d time_avg=%lf time_stddev=%lf num_runs=%d\n",
               n_all, np, omp_get_num_procs(), omp_get_max_threads(), average, variance, ITERATIONS_PER_SIZE);
        } else {
          printf("MIX Odd-Even MergeSort: N=%d np=%d num_proc=%d max_th=%d time_avg=%lf time_stddev=%lf num_runs=%d\n",
               n_all, np, omp_get_num_procs(), omp_get_max_threads(), average, variance, ITERATIONS_PER_SIZE);
        }
      }
      free(recv);
      free(buff);
      free(test);
    }
    free(order.data);
  }
  /////////////////////////////////////////////////////////////
  ////  E N D     O F     O R D E R I N G     C H O I C E  ////
  /////////////////////////////////////////////////////////////

  free(sizes);
  MPI_Finalize();
  return 0;
}