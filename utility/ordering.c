#include <stdlib.h>
#include <stdint.h>

struct MyVector {
  int32_t size;
  int32_t capacity;
  short* data;
};

struct MyVector GetVector(int32_t capacity) {
  if (capacity <= 0) {
    struct MyVector vec = {0, 0, NULL};
    return vec;
  }
  struct MyVector vec = {0, capacity, (short*)malloc(capacity * sizeof(short))};
  return vec;
}

void PushBack(struct MyVector* pvec, short value) {
  if (pvec->data == NULL) {
    ++pvec->capacity;
    pvec->data = (short*)malloc(sizeof(short));
  } else if (pvec->capacity == pvec->size) {
    pvec->capacity = (pvec->capacity * 3 + 1) >> 1;
    pvec->data = (short*)realloc(pvec->data, pvec->capacity * sizeof(short));
  }
  pvec->data[pvec->size] = value;
  ++pvec->size;
}


void SubCalcEvenOdd(struct MyVector** ord, short n) {
  ord[n] = (struct MyVector*)malloc(n * sizeof(struct MyVector));
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

struct MyVector* EvenOddOrdering(short n) {
  if (n < 1) {
    return NULL;
  } else if (n == 1) {
    struct MyVector* p_order = (struct MyVector*)malloc(sizeof(struct MyVector));
    *p_order = GetVector(0);
    return p_order;
  }
  struct MyVector** ordering = (struct MyVector**)malloc((n + 1) * sizeof(struct MyVector*));
  ordering[0] = NULL;

  ordering[1] = (struct MyVector*)malloc(sizeof(struct MyVector));
  ordering[1][0] = GetVector(0);

  ordering[2] = (struct MyVector*)malloc(2 * sizeof(struct MyVector));
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
  struct MyVector* result = ordering[n];
  free(ordering);
  return result;
}

struct MyVector* BitonicOrdering(short n) {
  if (n < 1) {
    return NULL;
  } else if (n == 1) {
    struct MyVector* p_order = (struct MyVector*)malloc(sizeof(struct MyVector));
    *p_order = GetVector(0);
    return p_order;
  }
  short p = 0, nn = 1;
  while (nn < n) {
    p += 1;
    nn <<= 1;
  }
  struct MyVector* ordering = (struct MyVector*)malloc(n * sizeof(struct MyVector));
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

void FreeOrdering(struct MyVector* ordering, short n) {
  for (short i = 0; i < n; ++i) {
    free(ordering[i].data);
  }
  free(ordering);
}