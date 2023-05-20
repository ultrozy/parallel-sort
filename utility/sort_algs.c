#include <stdint.h>

#define SWAP(type, temp, x, y) \
  type(temp) = (x);            \
  (x) = (y);                   \
  (y) = (temp)

typedef int ValueType;
typedef int32_t IndexType;

void HeapSort(ValueType *const array, const IndexType n) {
  for (IndexType i = (n >> 1) - 1; i != -1; --i) {
    IndexType index = i;

    // Heapify(array, array + n, index)
    while (index < n) {
      IndexType largest = index, left = (index << 1) + 1, right = left + 1;
      if (left < n && array[largest] < array[left]) {
        largest = left;
      }
      if (right < n && array[largest] < array[right]) {
        largest = right;
      }
      if (largest == index) {
        break;
      }
      SWAP(ValueType, temp, array[index], array[largest]);
      index = largest;
    }
  }
  for (IndexType i = n - 1; i != -1; --i) {
    SWAP(ValueType, temp, *array, array[i]);
    IndexType index = 0;

    // Heapify(array, array + i, 0)
    while (index < i) {
      IndexType largest = index, left = (index << 1) + 1, right = left + 1;
      if (left < i && array[largest] < array[left]) {
        largest = left;
      }
      if (right < i && array[largest] < array[right]) {
        largest = right;
      }
      if (largest == index) {
        break;
      }
      SWAP(, temp, array[index], array[largest]);
      index = largest;
    }
  }
}

ValueType *MergeSort(ValueType *array, ValueType *auxil, const IndexType n) {
  for (IndexType sort_len = 2; sort_len < (n << 1); sort_len <<= 1) {
    for (IndexType start1 = 0, start2 = sort_len >> 1; start1 < n; start1 += sort_len, start2 += sort_len) {
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

void QuickSort(ValueType *const array, const IndexType left, const IndexType right) {
  if (left + 1 >= right) {
    return;
  }
  ValueType pivot = array[(left + right) >> 1];
  IndexType li = left, ri = right - 1;
  IndexType le = li, re = ri;
  while (li <= ri) {
    while (array[li] < pivot) {
      ++li;
    }
    while (array[ri] > pivot) {
      --ri;
    }
    if (li >= ri) {
      break;
    }
    SWAP(ValueType, temp, array[li], array[ri]);
    if (array[li] == pivot) {
      SWAP(, temp, array[le], array[li]);
      ++le;
    }
    ++li;
    if (array[ri] == pivot) {
      SWAP(, temp, array[re], array[ri]);
      --re;
    }
    --ri;
  }

  if (li == ri) {
    ++li;
    --ri;
  }
  for (int k = left; k < le; ++k, --ri) {
    SWAP(ValueType, temp, array[k], array[ri]);
  }
  for (int k = right - 1; k > re; --k, ++li) {
    SWAP(ValueType, temp, array[k], array[li]);
  }
  QuickSort(array, left, ri + 1);
  QuickSort(array, li, right);
}

double MySqrt(double x) {
  double l = 1.0, r = x;
  if (l > r) {
    SWAP(double, temp, l, r);
  }
  while (r - l > 1e-4) {
    double m = 0.5 * (r + l), mm = m * m;
    if (mm <= x) {
      l = m;
    }
    if (mm >= x) {
      r = m;
    }
  }
  return l;
}