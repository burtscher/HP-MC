#ifndef CC_HELPER_H
#define CC_HELPER_H

#include "ECLgraph.h"
#include <vector>
#include <cstring>
#include <omp.h>

using std::vector;

static void color_in_component_csr(ECLgraph g, const int v, int& sum, int x,
  int &iter, vector<int> &candidate, const int valid_length,
  vector<int> &colors, int &visited_color, int* color_temp, int* prune) {
  //x is the pivot vertex
  iter++;
  color_temp[colors[v]] = iter;
  visited_color = 1;
  sum = 1;
  int candid_pointer = 0;
  int u = candidate[v];
  for (int j = g.nindex[u]; j < g.nindex[u + 1]; j++) {
    const int n = g.nlist[j];
    if (n == x) {
      visited_color = valid_length;
      break;
    }
    if (!prune[n]) {
      while ((candid_pointer < valid_length - 1) && (candidate[candid_pointer] < n ||
        candidate[candid_pointer] == -1)) candid_pointer++;
      if (candidate[candid_pointer] == n) {
        if (color_temp[colors[candid_pointer]] != iter) {
          color_temp[colors[candid_pointer]] = iter;
          visited_color++;
        }
      }
    }
  }
}

static void component(ECLgraph g, const int v, int& sum, vector<int> &label, int x, int &iter, int* wl1, int* wl2,
  int** adj, int* available, const int thread_index, const int valid_length, vector<int> &colors,
  int &visited_color, int* color_temp)
{
  // level-by-level BFS
  // x variable is the vertex we are not supposed to be adjacent to
  iter++;
  label[v] = iter;
  wl1[0] = v;
  color_temp[colors[v]] = iter;
  visited_color = 1;
  int g_level = 0;
  int size1 = 1;
  sum = 1;
  int size2 = 0;
  int level = 0;
  int stamp = iter;
  do {
    level++;
    size2 = 0;
    for (int i = 0; i < size1; i++) {
      const int u = wl1[i];
      for (int j = g.nindex[u]; j < g.nindex[u + 1]; j++) {
        const int n = g.nlist[j];
        if (available[n] == x || available[n] == -1 ||
          adj[x][available[n] - thread_index] == 1) continue;
	if (label[n] == 0) {
          label[n] = iter;
          wl2[size2++] = n;

	  if (color_temp[colors[n]] != iter) {
	    color_temp[colors[n]] = iter;
            visited_color++;
	  }
        }
      }
    }
    sum += size2;
    size1 = size2;
    if (size1 == 0) break;
    std::swap(wl1, wl2);
  } while (true);
}

static std::vector<int> component_count_csr(const ECLgraph g, int pivot_v, vector<int> &candidate,
  const int valid_length, vector<int> &colors,
   int* color_temp, int* prune, int &pivot_upper, int color_array_size) {
  vector<int> cc_color_size;
  int iter = 0;
  //std::fill(label.begin(), label.end(), 0);
  std::memset(color_temp, 0, sizeof(int) * (color_array_size + 1));
  for (int i = 0; i < valid_length; i++) {
    int v = candidate[i];
    if (v == -1) {
      cc_color_size.push_back(0);
      continue;
    }
      //continue;
    int sum;
    int visited_colors;
    color_in_component_csr(g, i, sum, pivot_v, iter,
      candidate, valid_length, colors, visited_colors,
      color_temp, prune);
    cc_color_size.push_back(visited_colors);
    if (v == pivot_v) pivot_upper = visited_colors;
  }
  return cc_color_size;
}

static std::vector<int> component_count(const ECLgraph g, int pivot_v, vector<int> &label,
  int* wl1, int* wl2, int** adj, int* available, const int thread_index,
  const int valid_length, vector<int> &colors, vector<int> &cc_color_size, int* color_temp) {
  int iter = 0;
  std::fill(label.begin(), label.begin() + g.nodes, 0);
  std::memset(color_temp, 0, sizeof(int) * g.nodes);
  cc_color_size.clear();
  std::vector<int> cc_size;
  for (int i = 0; i < g.nodes; i++) {
    int sum = 0;
    if (available[i] == pivot_v || available[i] == -1 || adj[available[i]][pivot_v - thread_index] == 1) continue;
    if (label[i] == 0) {
      int visited_colors;
      component(g, i, sum, label, pivot_v, iter, wl1,
        wl2, adj, available, thread_index, valid_length, colors, visited_colors, color_temp);
      cc_size.push_back(sum);
      cc_color_size.push_back(visited_colors);
    }
  }
  return cc_size;
}

#endif //CC_HELPER_H

