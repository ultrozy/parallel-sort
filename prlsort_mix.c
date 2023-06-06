#include <mpi.h>
#include <omp.h>
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
  1) MyVector               - std::vector analogue with limited functionality
  2) SortNet                - Array of MyVectors. Technically contains sorting network
  3) GetVector              - returns MyVector constructed by default
  4) OEmergenet, OEsortnet  - subroutines for OddEvenOrdering
  5) OddEvenOrdering        - returns Odd-Even ordering
  6) BitonicOrdering        - returns Bitonic  ordering
*/
struct MyVector {
  int32_t size;
  int32_t capacity;
  short *data;
};
typedef struct MyVector* SortNet;

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

void OEmergenet(SortNet* all_mergenets, short n) {
  all_mergenets[n] = (SortNet)malloc(n * sizeof(struct MyVector));
  SortNet net = all_mergenets[n];
  for (short i = 0; i < n; ++i) {
    net[i] = GetVector(0);
  }
  short a_size = ((n + 1) >> 1), b_size = n - a_size;
  short a_odd = (a_size >> 1), b_odd = (b_size >> 1), odds = a_odd + b_odd;
  short a_even = a_size - a_odd, b_even = b_size - b_odd, evens = a_even + b_even;

  // MERGE EVEN
  //// A
  for (short i = 0; i < a_size; i += 2) {
    short converted = i >> 1;
    for (int32_t j = 0; j < all_mergenets[evens][converted].size; ++j) {
      short elem = all_mergenets[evens][converted].data[j];
      PushBack(&net[i], elem < a_even ? elem << 1 : a_size + ((elem - a_even) << 1));
    }
  }
  //// B
  for (short i = a_size; i < n; i += 2) {
    short converted = ((i - a_size) >> 1) + a_even;
    for (int32_t j = 0; j < all_mergenets[evens][converted].size; ++j) {
      short elem = all_mergenets[evens][converted].data[j];
      PushBack(&net[i], elem < a_even ? elem << 1 : a_size + ((elem - a_even) << 1));
    }
  }
  // MERGE ODD
  //// A
  for (short i = 1; i < a_size; i += 2) {
    short converted = i >> 1;
    for (int32_t j = 0; j < all_mergenets[odds][converted].size; ++j) {
      short elem = all_mergenets[odds][converted].data[j];
      PushBack(&net[i], elem < a_odd ? (elem << 1) + 1 : a_size + 1 + ((elem - a_odd) << 1));
    }
  }
  //// B
  for (short i = a_size + 1; i < n; i += 2) {
    short converted = ((i - a_size) >> 1) + a_odd;
    for (int32_t j = 0; j < all_mergenets[odds][converted].size; ++j) {
      short elem = all_mergenets[odds][converted].data[j];
      PushBack(&net[i], elem < a_odd ? (elem << 1) + 1 : a_size + 1 + ((elem - a_odd) << 1));
    }
  }

  // EMERGE
  for (short i = 1; i < n - 1; i += 2) {
    PushBack(&net[i], i + 1);
    PushBack(&net[i + 1], i);
  }
}
void OEsortnet(SortNet* all_sortnets, SortNet* all_mergenets, short n) {
  all_sortnets[n] = (SortNet)malloc(n * sizeof(struct MyVector));
  short a_size = ((n + 1) >> 1), b_size = n - a_size;
  SortNet net = all_sortnets[n], merge = all_mergenets[n];
  SortNet asort = all_sortnets[a_size], bsort = all_sortnets[b_size];

  // A
  for (short i = 0; i < a_size; ++i) {
    net[i] = GetVector(asort[i].size + merge[i].size);
    for (int32_t j = 0; j < asort[i].size; ++j) {
      PushBack(&net[i], asort[i].data[j]);
    }
    for (int32_t j = 0; j < merge[i].size; ++j) {
      PushBack(&net[i], merge[i].data[j]);
    }
  }

  // B
  for (short i = a_size; i < n; ++i) {
    net[i] = GetVector(bsort[i - a_size].size + merge[i].size);
    for (int32_t j = 0; j < bsort[i - a_size].size; ++j) {
      PushBack(&net[i], bsort[i - a_size].data[j] + a_size);
    }
    for (int32_t j = 0; j < merge[i].size; ++j) {
      PushBack(&net[i], merge[i].data[j]);
    }
  }
}
SortNet OddEvenOrdering(short n) {
  if (n < 1) {
    return NULL;
  } else if (n == 1) {
    SortNet sortnet = (SortNet)malloc(sizeof(struct MyVector));
    sortnet[0] = GetVector(0);
    return sortnet;
  }

  // Sort(k)
  SortNet* all_sortnets = (SortNet*)malloc((n + 1) * sizeof(SortNet));
  all_sortnets[0] = NULL;
  all_sortnets[1] = (SortNet)malloc(sizeof(struct MyVector));
  all_sortnets[1][0] = GetVector(0);

  // [k] -> Merge(ceil(k/2), floor(k/2))
  SortNet* all_mergenets = (SortNet*)malloc((n + 1) * sizeof(SortNet));
  all_mergenets[0] = NULL;
  all_mergenets[1] = (SortNet)malloc(sizeof(struct MyVector));
  all_mergenets[1][0] = GetVector(0);
  all_mergenets[2] = (SortNet)malloc(2 * sizeof(struct MyVector));
  all_mergenets[2][0] = GetVector(1);
  all_mergenets[2][1] = GetVector(1);
  PushBack(&all_mergenets[2][0], 1);
  PushBack(&all_mergenets[2][1], 0);

  for (short k = 3; k <= n; ++k) {
    OEmergenet(all_mergenets, k);
  }
  for (short k = 2; k <= n; ++k) {
    OEsortnet(all_sortnets, all_mergenets, k);
  }

  SortNet result = all_sortnets[n];
  for (short k = 1; k < n; ++k) {
    for (short j = 0; j < k; ++j) {
      free(all_sortnets[k][j].data);
      free(all_mergenets[k][j].data);
    }
    free(all_sortnets[k]);
    free(all_mergenets[k]);
  }
  for (short j = 0; j < n; ++j) {
    free(all_mergenets[n][j].data);
  }
  free(all_mergenets[n]);
  free(all_mergenets);
  free(all_sortnets);
  return result;
}
SortNet BitonicOrdering(short n) {
  if (n < 1) {
    return NULL;
  } else if (n == 1) {
    SortNet sortnet = (SortNet)malloc(sizeof(struct MyVector));
    sortnet[0] = GetVector(0);
    return sortnet;
  }
  short p = 0, nn = 1;
  while (nn < n) {
    p += 1;
    nn <<= 1;
  }
  SortNet sortnet = (SortNet)malloc(n * sizeof(struct MyVector));
  for (short i = 0; i < n; ++i) {
    sortnet[i] = GetVector((p * (p - 1)) >> 1);
  }
  for (short halflen = 1, len = 2; halflen < nn; halflen <<= 1, len <<= 1) {
    // First Block
    for (short start = 0; start < n; start += len) {
      for (short i = 0; i < halflen; ++i) {
        short j = len - i - 1;
        if (i + start < n && j + start < n) {
          PushBack(&sortnet[start + i], start + j);
          PushBack(&sortnet[start + j], start + i);
        }
      }
    }
    // Other Blocks
    for (short dist = halflen >> 1, delta = halflen; dist > 0; dist >>= 1, delta >>= 1) {
      for (short start = 0; start < n; start += delta) {
        for (short i = start; i < start + dist; ++i) {
          if (i + dist < n) {
            PushBack(&sortnet[i], i + dist);
            PushBack(&sortnet[i + dist], i);
          }
        }
      }
    }
  }
  return sortnet;
}

ValueType *MergeSortOMP(ValueType *array, ValueType *auxil, const IndexType size) {
  // SORTED LENGTH LOOP. sequence of sort_len will come up AFTER the iteration
  for (IndexType sort_len = 2; sort_len < (size << 1); sort_len <<= 1) {
    // MERGE LOOP
#pragma omp parallel for
    for (IndexType start1 = 0; start1 < size; start1 += sort_len) {
      IndexType start2 = start1 + (sort_len >> 1);
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

// The only difference with MergeSortMPI: instead of MergeSort - MergeSortOMP is used
void MergeSortMIX() {
  // INPUT
  IndexType size, capacity;
  ValueType *array;
  InputForMPI(array, &size, &capacity);

  // From now on capacity is constant

  // SORT own part and fill VALUETYPE_MAX
  ValueType *auxil = (ValueType *)malloc(capacity * sizeof(ValueType));
  if (MergeSortOMP(array, auxil, size) != array) {
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
  struct MyVector *ordering = OddEvenOrdering(np);
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
