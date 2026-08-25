#ifndef INTERSECT_H
#define INTERSECT_H
#include <climits>
static inline void common_matrix(const int* const __restrict__ list1, const int end1, const int* const __restrict__ list2,
  const int end2, int** const __restrict__ Adj, int cur_length, int* const __restrict__ count_array, const int thread_index)
{
  //update Adj
  int val1, val2 = (end2 > 0) ? list2[0] : INT_MAX;  // list2 may be empty
  int pos1 = 0, pos2 = 0;
  while (true) {
    while ((pos1 < end1) && ((val1 = list1[pos1]) < val2)) {
      pos1++;
    }
    if (pos1 == end1) break;
    if (val2 < val1) {
      pos2++;
      while ((pos2 < end2) && ((val2 = list2[pos2]) < val1)) {
        pos2++;
      }
      if (pos2 == end2) break;
    }
    if (val1 == val2) {
      //#pragma omp atomic write
      Adj[thread_index + cur_length][pos1] = 1;
      //#pragma omp atomic write
      Adj[thread_index + pos1][cur_length] = 1;
      //#pragma omp atomic update
      count_array[thread_index + cur_length]++;
      //#pragma omp atomic update
      count_array[thread_index + pos1]++;

      pos2++;
      if (pos2 == end2) break;
      val2 = list2[pos2];
    }
    pos1++;
  }
}

#endif //INTERSECT_H

