#include <stdlib.h>
#include <stdint.h>

#define SWAP(type, temp, x, y) \
  type(temp) = (x);            \
  (x) = (y);                   \
  (y) = (temp)

typedef int ValueType;
typedef int32_t IndexType;

int *GetSortedInts(const IndexType n) {
  int *array = (int *)malloc(n * sizeof(int));
  for (IndexType i = 0; i < n; ++i) {
    array[i] = i - (n >> 1);
  }
  return array;
}

int *GetRandomInts(const IndexType n, const IndexType min, const IndexType max, unsigned int seed) {
  if (max < min) {
    return NULL;
  }
  srand(seed);
  int *array = (int *)malloc(n * sizeof(int));
  for (IndexType i = 0; i < n; ++i) {
    array[i] = min + rand() % (max - min + 1);
  }
  return array;
}

void SetRandomInts(int* array, const IndexType n, const IndexType min, const IndexType max, unsigned int seed) {
  if (max < min) {
    return;
  }
  srand(seed);
  for (IndexType i = 0; i < n; ++i) {
    array[i] = min + rand() % (max - min + 1);
  }
}

ValueType *GetReverseTest(ValueType *const array, const IndexType n) {
  for (IndexType i = 0; i < (n >> 1); ++i) {
    SWAP(ValueType, temp, array[i], array[n - 1 - i]);
  }
  return array;
}

ValueType *GetRandomizeTest(ValueType *const array, const IndexType n, unsigned int seed) {
  if (n <= 0) {
    return array;
  }
  srand(seed);
  for (IndexType i = n - 1; i > 0; --i) {
    IndexType j = rand() % (i + 1);
    if (j < i) {
      SWAP(ValueType, temp, array[i], array[j]);
    }
  }
  return array;
}

ValueType *GetSlightlyRandomizedTest(ValueType *array, const IndexType n, IndexType window, unsigned int seed) {
  if (n <= 0) {
    return array;
  }
  if (window <= 0) {
    window = 1;
  }
  IndexType *indices = (IndexType *)malloc(n * sizeof(IndexType));
  for (IndexType i = 0; i < n; ++i) {
    indices[i] = i;
  }

  srand(seed);
  for (IndexType i = n - 1; i > 0; --i) {
    IndexType start = (indices[i] < window ? 0 : indices[i] - window);
    IndexType j = start + rand() % (i + 1 - start);
    if (j < i) {
      SWAP(IndexType, temp, indices[i], indices[j]);
    }
  }

  ValueType *new_array = (ValueType *)malloc(n * sizeof(ValueType));
  for (IndexType i = 0; i < n; ++i) {
    new_array[indices[i]] = array[i];
  }
  free(array);
  free(indices);
  return new_array;
}
