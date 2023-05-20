#include <mpi.h>
#include <stdlib.h>

#define SWAP(type, temp, x, y) type(temp) = (x); (x) = (y); (y) = (temp)

typedef int ValueType;
typedef int IndexType;
#define VALUETYPE_MAX INT32_MAX;
#define MPI_VALUETYPE MPI_INT
#define MPI_INDEXTYPE MPI_INT

#define USE_ODDEVEN_ORDERING
// #define USE_BITONIC_ORDERING

// For MPI to work input and output for single process must be specified
// Make sure, that CAPACITY IS EQUAL for all processes
void InputForMPI(ValueType *array, IndexType *psize, IndexType *pcapacity);
void OutputForMPI(ValueType *array, IndexType size);

/*
Taken from utility/ordering.c:
  1) MyVector         - std::vector analogue with limited functionality 
  2) GetVector        - returns MyVector constructed by default
  3) SubCalcEvenOdd   - subroutine for EvenOddOrdering
  4) EvenOddOrdering  - returns Odd-Even ordering
  5) BitonicOrdering  - returns Bitonic  ordering
*/
struct MyVector {
  int32_t size;
  int32_t capacity;
  short *data;
};

struct MyVector GetVector(int32_t capacity) {
  if (capacity <= 0) {
    struct MyVector vec = {0, 0, NULL};
    return vec;
  }
  struct MyVector vec = {0, capacity, (short *)malloc(capacity * sizeof(short))};
  return vec;
}
void PushBack(struct MyVector *pvec, short value) {
  if (pvec->data == NULL) {
    ++pvec->capacity;
    pvec->data = (short *)malloc(sizeof(short));
  } else if (pvec->capacity == pvec->size) {
    pvec->capacity = (pvec->capacity * 3 + 1) >> 1;
    pvec->data = (short *)realloc(pvec->data, pvec->capacity * sizeof(short));
  }
  pvec->data[pvec->size] = value;
  ++pvec->size;
}

void SubCalcEvenOdd(struct MyVector **ord, short n) {
  ord[n] = (struct MyVector *)malloc(n * sizeof(struct MyVector));
  for (short i = 0; i < n; ++i) {
    ord[n][i] = GetVector(0);
  }

  // FIRST PHASE
  // Merge A
  short a_size = ((n + 1) >> 1), b_size = n - a_size;
  for (short i = 0; i < a_size; ++i) {
    for (int32_t j = 0; j < ord[a_size][i].size; ++j) {
      PushBack(&ord[n][i], ord[a_size][i].data[j]);
    }
  }

  // Merge B
  for (short i = a_size; i < n; ++i) {
    for (int32_t j = 0; j < ord[b_size][i - a_size].size; ++j) {
      short elem = ord[b_size][i - a_size].data[j];
      PushBack(&ord[n][i], elem + a_size);
    }
  }

  // SECOND PHASE
  short a_odd = (a_size >> 1), b_odd = (b_size >> 1), odds = a_odd + b_odd;
  short a_even = a_size - a_odd, b_even = b_size - b_odd, evens = a_even + b_even;
  // Evens
  //// A
  for (short i = 0; i < a_size; i += 2) {
    short converted = i >> 1;
    for (int32_t j = 0; j < ord[evens][converted].size; ++j) {
      short elem = ord[evens][converted].data[j];
      PushBack(&ord[n][i], elem < a_even ? elem << 1 : a_size + ((elem - a_even) << 1));
    }
  }
  //// B
  for (short i = a_size; i < n; i += 2) {
    short converted = ((i - a_size) >> 1) + a_even;
    for (int32_t j = 0; j < ord[evens][converted].size; ++j) {
      short elem = ord[evens][converted].data[j];
      PushBack(&ord[n][i], elem < a_even ? elem << 1 : a_size + ((elem - a_even) << 1));
    }
  }
  // Odds
  //// A
  for (short i = 1; i < a_size; i += 2) {
    short converted = i >> 1;
    for (int32_t j = 0; j < ord[odds][converted].size; ++j) {
      short elem = ord[odds][converted].data[j];
      PushBack(&ord[n][i], elem < a_odd ? (elem << 1) + 1 : a_size + 1 + ((elem - a_odd) << 1));
    }
  }
  //// B
  for (short i = a_size + 1; i < n; i += 2) {
    short converted = ((i - a_size) >> 1) + a_odd;
    for (int32_t j = 0; j < ord[odds][converted].size; ++j) {
      short elem = ord[odds][converted].data[j];
      PushBack(&ord[n][i], elem < a_odd ? (elem << 1) + 1 : a_size + 1 + ((elem - a_odd) << 1));
    }
  }

  // THIRD PHASE
  for (int i = 1; i < n - 1; i += 2) {
    PushBack(&ord[n][i], i + 1);
    PushBack(&ord[n][i + 1], i);
  }
}
struct MyVector *EvenOddOrdering(short n) {
  if (n < 1) {
    return NULL;
  } else if (n == 1) {
    struct MyVector *p_order = (struct MyVector *)malloc(sizeof(struct MyVector));
    *p_order = GetVector(0);
    return p_order;
  }
  struct MyVector **ordering = (struct MyVector **)malloc((n + 1) * sizeof(struct MyVector *));
  ordering[0] = NULL;

  ordering[1] = (struct MyVector *)malloc(sizeof(struct MyVector));
  ordering[1][0] = GetVector(0);

  ordering[2] = (struct MyVector *)malloc(2 * sizeof(struct MyVector));
  ordering[2][0] = GetVector(1);
  ordering[2][1] = GetVector(1);
  PushBack(&ordering[2][0], 1);
  PushBack(&ordering[2][1], 0);

  for (short i = 3; i <= n; ++i) {
    SubCalcEvenOdd(ordering, i);
  }
  for (short i = 1; i < n; ++i) {
    for (short j = 0; j < i; ++j) {
      free(ordering[i][j].data);
    }
    free(ordering[i]);
  }
  struct MyVector *result = ordering[n];
  free(ordering);
  return result;
}
struct MyVector *BitonicOrdering(short n) {
  if (n < 1) {
    return NULL;
  } else if (n == 1) {
    struct MyVector *p_order = (struct MyVector *)malloc(sizeof(struct MyVector));
    *p_order = GetVector(0);
    return p_order;
  }
  short p = 0, nn = 1;
  while (nn < n) {
    p += 1;
    nn <<= 1;
  }
  struct MyVector *ordering = (struct MyVector *)malloc(n * sizeof(struct MyVector));
  for (short i = 0; i < n; ++i) {
    ordering[i] = GetVector((p * (p - 1)) >> 1);
  }
  for (short halflen = 1, len = 2; halflen < nn; halflen <<= 1, len <<= 1) {
    // First Block
    for (short start = 0; start < n; start += len) {
      for (short i = 0; i < halflen; ++i) {
        short j = len - i - 1;
        if (i + start < n && j + start < n) {
          PushBack(&ordering[start + i], start + j);
          PushBack(&ordering[start + j], start + i);
        }
      }
    }
    // Other Blocks
    for (short dist = halflen >> 1, delta = halflen; dist > 0; dist >>= 1, delta >>= 1) {
      for (short start = 0; start < n; start += delta) {
        for (short i = start; i < start + dist; ++i) {
          if (i + dist < n) {
            PushBack(&ordering[i], i + dist);
            PushBack(&ordering[i + dist], i);
          }
        }
      }
    }
  }
  return ordering;
}

ValueType *MergeSort(ValueType *array, ValueType *auxil, const IndexType size) {
  // SORTED LENGTH LOOP. sequence of sort_len will come up AFTER the iteration
  for (IndexType sort_len = 2; sort_len < (size << 1); sort_len <<= 1) {
    // MERGE LOOP
    for (IndexType start1 = 0, start2 = sort_len >> 1; start1 < size; start1 += sort_len, start2 += sort_len) {
      IndexType end1 = start2, end2 = start1 + sort_len, i1 = start1, i2 = start2, j = start1;
      if (end1 > size) {
        end1 = size;
      }
      if (end2 > size) {
        end2 = size;
      }
      // MERGE
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
      // FINISH MERGE
      while (i1 < end1) {
        auxil[j++] = array[i1++];
      }
      while (i2 < end2) {
        auxil[j++] = array[i2++];
      }
    }
    // To avoid memcpy between (array) and (auxil) SWAP pointers
    SWAP(ValueType *, temp, array, auxil);
  }
  return array;
}

// The only difference with MergeSortMIX: instead of MergeSortOMP - MergeSort is used
void MergeSortMPI() {
  // INPUT
  IndexType size, capacity;
  ValueType *array;
  InputForMPI(array, &size, &capacity);

  // From now on capacity is constant

  // SORT own part and fill VALUETYPE_MAX
  ValueType *auxil = (ValueType *)malloc(capacity * sizeof(ValueType));
  if (MergeSort(array, auxil, size) != array) {
    SWAP(ValueType *, temp, array, auxil);
  }
  for (IndexType i = size; i < capacity; ++i) {
    array[i] = VALUETYPE_MAX;
  }

  MPI_Init(NULL, NULL);
  int id, np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

// DEFINE ORDERING
#ifdef USE_ODDEVEN_ORDERING
  struct MyVector *ordering = EvenOddOrdering(np);
#else
  struct MyVector *ordering = BitonicOrdering(np);
#endif
  struct MyVector order = ordering[id];
  for (short i = 0; i < np; ++i) {
    if (i != id) {
      free(ordering[i].data);
    }
  }
  free(ordering);

  // MERGE
  ValueType *recvb = (ValueType *)malloc(capacity * sizeof(ValueType));
  for (short iord = 0; iord < order.size; ++iord) {
    MPI_Sendrecv( array, capacity, MPI_VALUETYPE, order.data[iord], 0,
                  recvb, capacity, MPI_VALUETYPE, order.data[iord], 0,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (id < order.data[iord]) {
      IndexType i = 0, j = 0;
      for (IndexType k = 0; k < capacity; ++k) {
        if (recvb[j] < array[i]) {
          auxil[k] = recvb[j++];
        } else {
          auxil[k] = array[i++];
        }
      }
    } else {
      IndexType i = capacity - 1, j = capacity - 1;
      for (IndexType k = capacity - 1; k != -1; --k) {
        if (array[i] < recvb[j]) {
          auxil[k] = recvb[j--];
        } else {
          auxil[k] = array[i--];
        }
      }
    }
    SWAP(ValueType *, temp, array, auxil);
  }
  free(recvb);
  free(auxil);

  // Gathering full array size into (num_output)
  IndexType num_output;
  MPI_Allreduce(&size, &num_output, 1, MPI_INDEXTYPE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();

  // Determine how many elements to write
  num_output -= id * capacity;
  if (num_output > capacity) {
    num_output = capacity;
  } else if (num_output < 0) {
    num_output = 0;
  }

  OutputForMPI(array, num_output);
  free(array);
}
