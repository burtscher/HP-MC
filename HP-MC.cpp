/*
HP-MC is a C++/OpenMP code for quickly finding the exact maximum clique of large sparse graphs.

Copyright (c) 2026, Cameron Bradley, and Martin Burtscher

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

URL: The latest version of this code is available at https://github.com/CameronBradley1/HP-MC

Publication: This work is described in detail in the following paper.


Sponsor: This work has been supported by the National Science Foundation under Award #1955367.
*/


#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cassert>
#include <sys/time.h>
#include "ECLgraph.h"
#include "ECLgc.h"
#include "cc_helper.h"
#include <parallel/algorithm>
#include <omp.h>
#include <iostream>
#include <climits>
#include <cstring>
#include <cmath>
#include "intersect.h"
#include <execution>
#include <unordered_set>
#include <sys/mman.h>

/*
 The input graph must meet the following conditions, which are not checked:
 - The graph cannot contain self-edges
 - The graph cannot contain multiple edges between the same pair of vertices
 - The graph must be undirected, i.e., every edge must appear twice in the CSR format, once in each direction
 - The individual adjacency lists must be sorted in increasing order
*/
using std::vector;
static int num_threads;
static int adj_size = 0;
int maxed;
static int* index_map;
static int* clique_embed;
static int max_size = 0;
int max_wl;
int idle_threads = 0;
int*** g_available;
int*** g_r_ed;
int** g_wl1;
int** g_wl2;
int** g_color_temp;
int* count_array;
int* prune;
int** word_Adj;
bool* task_ready;
int* thread_depth;
ECLgraph* global_g_temp;
static const int MAXN = 2048;
using U64 = uint64_t;
static const int W = (MAXN + 63) / 64;
struct Bitset {
  U64 a[W];

  Bitset() { memset(a, 0, sizeof(a)); }

  inline void set(int i) { a[i >> 6] |= (1ULL << (i & 63)); }
  inline void reset(int i) { a[i >> 6] &= ~(1ULL << (i & 63)); }
  inline bool test(int i) const { return a[i >> 6] & (1ULL << (i & 63)); }

  inline void operator&=(const Bitset& o) {
    for (int i = 0; i < W; i++) a[i] &= o.a[i];
  }

  inline void operator^=(const Bitset& o) {
    for (int i = 0; i < W; i++) a[i] ^= o.a[i];
  }

  inline int popcount() const {
    int s = 0;
    for (int i = 0; i < W; i++) s += __builtin_popcountll(a[i]);
    return s;
  }

  inline bool empty() const {
    for (int i = 0; i < W; i++) if (a[i]) return false;
    return true;
  }

  inline int pick() const {
    for (int i = 0; i < W; i++) {
      if (a[i]) return (i << 6) + __builtin_ctzll(a[i]);
    }
    return -1;
  }

};

inline int get_indices(const Bitset& bs, vector<int> &out, int tid,
vector<int> &color, const int stamp) {
  out.clear();
  int cnt = 0;
  //Bitset color_temp;
  for (int i = 0; i < W; i++) {
    U64 x = bs.a[i];
    while (x) {
      int b = __builtin_ctzll(x);
      int u = (i << 6) + b;
      if (g_color_temp[tid][color[u]] != stamp) cnt++;
      g_color_temp[tid][color[u]] = stamp;
      //*/
      out.push_back(u);
      x &= (x - 1);
    }
  }
  return cnt;
}
struct Global_task
{
  int v;
  int degree;
  int effective;
  int major = 0;
  int lower_clique = 0;
  vector<int> partial_clique;

  Global_task(const int vertex, const int degree, const int effective) :
  v(vertex), degree(degree), effective(effective) {}
};


struct CPUTimer
{
  timeval beg, end;
  CPUTimer() {}
  ~CPUTimer() {}
  void start() {gettimeofday(&beg, NULL);}
  double elapsed() {gettimeofday(&end, NULL); return end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) / 1000000.0;}
};

CPUTimer timer;

vector<vector<vector<int>>> color(64,
  vector<vector<int>>(MAXN,
    vector<int>(MAXN)
  )
);

int matrix_BRB(const ECLgraph g, int* const __restrict__ ed, int start_size,
  int latest_level, int** const __restrict__ adj, int** available, int valid_length, int** r_ed,
  int& current_max, const int thread_index, int tid, ECLgraph g_temp,
  int* color_temp, int* wl1, int* wl2, bool task_allowed, const int stolen_index, vector<int> &colors) {
  int upper_bound = start_size + valid_length;
  int new_level = latest_level + 1;
  int pivot = -1;

  if (max_size + 1 == maxed + 1) {
    return max_size;
  }

  thread_depth[tid] = new_level;

  bool change;
  int length = 0;

  if (start_size > max_size) {
    upper_bound = start_size;
    current_max = start_size;
      max_size = start_size;
  }

  if (valid_length == 0) {
    thread_depth[tid] = latest_level;
    return start_size;
  }

  int total_cols = -1;
  g_temp.nodes = valid_length;
  g_temp.nindex[0] = 0;
  int edge_number = 0;
  for (int x = 0; x < valid_length; x++) {
    const int v_x = available[latest_level][x];
    for (int y = 0; y < valid_length; y++) {
      const int v_y = available[latest_level][y];
      if (x == y) continue;
      if (adj[v_x][v_y - thread_index] == 1) {
        g_temp.nlist[edge_number++] = y;
      }
    }
    g_temp.nindex[x + 1] = edge_number;
  }
  g_temp.edges = edge_number;

  vector<int> idx(valid_length);
  for (int i = 0; i < valid_length; i++) {
    idx[i] = i;
  }

  memset(color_temp, 0, sizeof(int) * (valid_length + 1));
  int stamp = -1;
  length = 0;
  int pivot_ub;
  int piv_deg;
  int high_color_v;
  do {
    high_color_v = 0;
    piv_deg = 0;
    change = false;
    int max_color = 1;
    stamp = 1;
    pivot_ub = 0;
    for (int x = 0; x < valid_length; x++) {
      stamp++;
      if (available[latest_level][x] == -1) continue;
      const int check = available[latest_level][x];
      if (prune[index_map[check]]) {
        available[latest_level][x] = -1;
      }
      for (int y = g_temp.nindex[x]; y < g_temp.nindex[x + 1]; y++) {
        int v = g_temp.nlist[y];
        if (v >= x) break;
        if (available[latest_level][v] == -1) continue;
        color_temp[colors[v]] = stamp;
      }
      int c = 1;
      while (color_temp[c] == stamp) {
        ++c;
      }
      colors[x] = c;
      high_color_v += (c >= max_size - start_size);
      max_color = std::max(max_color,c);
    }
    total_cols = max_color;

    if (total_cols + start_size <= max_size) {
      thread_depth[tid] = latest_level;
      return total_cols + start_size;
    }

    upper_bound = start_size + total_cols;
    length = 0;
    stamp = -1;
    piv_deg = -1;
    pivot = -1;
    for (int i = 0; i < valid_length; i++) {
      if (available[latest_level][i] == -1) continue;
      int v = i;
      int cnt = 0;
      int deg = 0;
      const int v_beg = g_temp.nindex[v];
      const int v_end = g_temp.nindex[v + 1];
      for (int j = v_beg; j < v_end; j++) {
        const int v_y = g_temp.nlist[j];
        const int og_vy = available[latest_level][v_y];
        if (available[latest_level][v_y] == -1) continue;
        if (color_temp[colors[v_y]] != stamp && !prune[index_map[og_vy]]) {
          color_temp[colors[v_y]] = stamp;
          cnt++;
        }
        deg += (!prune[index_map[og_vy]]);
      }
      if (cnt + start_size >= max_size && !prune[index_map[available[latest_level][i]]]) {
        if (deg > piv_deg) {
	        pivot = idx[i]; piv_deg = deg; pivot_ub = cnt;
	      }
	      length++;
      } else {
        available[latest_level][v] = -1;
        change = true;
      }
      stamp--;
    }
  } while (change && length + start_size > max_size);

  if (total_cols + start_size <= max_size || length + start_size <= max_size) {
    thread_depth[tid] = latest_level;
    return upper_bound;
  }

  if (pivot == -1) {
    thread_depth[tid] = latest_level;
    return upper_bound;
  }

  length = 0;
  for (int i = 0; i < valid_length; i++) {
    if (available[latest_level][i] != -1) {
      idx[length++] = idx[i];
    }
  }
  valid_length = length;

  const int pivot_vertex = available[latest_level][pivot];
  std::vector<int> skipped_v(valid_length, 0);
  vector<bool> removed_cc_labels(valid_length, false);
  vector<int> cc_remaining(valid_length, 0);
  vector<int>label_v(g_temp.nodes);
  const int start_const = start_size;
  vector<int> cc_color_size;
  std::vector<int> cc_size = component_count(g_temp, pivot_vertex, label_v, wl1, wl2, adj, available[latest_level],
    thread_index, valid_length, colors, cc_color_size, color_temp);

  vector<int> next_color(valid_length);
  int pivot_ub2 = 0;
  stamp = INT_MAX;
  int next_length = 0;
  for (int i = 0; i < valid_length; i++) {
    if (idx[i] == pivot) continue;
    int other_v = available[latest_level][idx[i]];
    if (other_v != -1 && adj[other_v][pivot_vertex - thread_index] == 1 &&
    !prune[index_map[other_v]]
    ) {
      if (color_temp[colors[idx[i]]] != stamp) {
        color_temp[colors[idx[i]]] = stamp;
        pivot_ub2++;
      }
      next_color[next_length] = colors[idx[i]];
      available[new_level][next_length++] = other_v;
    }
  }

  int potential_upper = pivot_ub2 + start_size;


  if (potential_upper >= max_size) {
    upper_bound = matrix_BRB(g, ed, start_const + 1, new_level, adj, available, next_length, r_ed,
              max_size, thread_index, tid, g_temp, color_temp, wl1, wl2,
              task_allowed, stolen_index, next_color);
    thread_depth[tid] = new_level;
  }

  for (int i = 0; i < cc_size.size(); i++) {
    if (cc_color_size[i] == 1) {
      removed_cc_labels[i + 1] = true;
      continue;
    }
    if (cc_color_size[i] + (potential_upper) <= max_size) {
      removed_cc_labels[i + 1] = true;
    }
    else {
      cc_remaining[i] = 1;
    }
  }

  int saved_max = max_size;

  for (int i = 0; i < valid_length; i++) {
    int vertex;
    vertex = available[latest_level][idx[i]];
    if (vertex == -1 || prune[index_map[vertex]]) continue;
    if ((vertex == pivot_vertex)) continue;
    if ((adj[pivot_vertex][vertex - thread_index] == 1)) continue;
    if (removed_cc_labels[label_v[idx[i]]]) {
      available[latest_level][idx[i]] = -1;
      if (latest_level == 0) prune[index_map[vertex]] = true;
      continue;
    }
    if (cc_remaining[label_v[idx[i]] - 1] > 0) {
      skipped_v[i] = 1;
      cc_remaining[label_v[idx[i]] - 1]--;
      continue;
    }

    #pragma omp task if((valid_length >= adj_size && idle_threads > 0)) priority(adj_size + latest_level) \
    firstprivate(thread_index, i, vertex, latest_level, saved_max) \
    shared(g_color_temp, colors, idle_threads, upper_bound, global_g_temp, idx, valid_length, start_size, \
    start_const, tid, pivot_vertex, new_level, skipped_v, task_allowed, g_available, g_r_ed)
    {
      int stealer = omp_get_thread_num();
      if (stealer != tid) {
        #pragma omp atomic update
        idle_threads--;
      }
      memset(g_color_temp[stealer], 0, sizeof(int) * (valid_length + 1));
      vector<int> task_color(valid_length);
      int steal_index = stealer * adj_size;
      int task_length = 0;
      int stolen_max = saved_max;
      int task_ub = 0;
      int task_stamp = 1;
      int saved_depth = thread_depth[stealer];
      for (int j = 0; j < valid_length; j++) {
        int other_v = g_available[tid][latest_level][idx[j]];
        if (other_v == -1 || j == i || removed_cc_labels[label_v[idx[j]]]) continue;
        if (adj[other_v][vertex - thread_index] == 1 && !prune[index_map[other_v]]) {
          task_ub += (g_color_temp[stealer][colors[idx[j]]] != task_stamp);
          g_color_temp[stealer][colors[idx[j]]] = task_stamp;
          task_color[task_length] = colors[idx[j]];
          g_available[stealer][saved_depth][task_length++] = other_v;
        }
      }
      int maximal = task_ub + start_const + 1;
      if (maximal > stolen_max) {
        maximal = matrix_BRB(g, ed, start_const + 1, saved_depth, adj,
          g_available[stealer], task_length, g_r_ed[stealer],max_size, thread_index, stealer,
          global_g_temp[stealer], g_color_temp[stealer], g_wl1[stealer], g_wl2[stealer], task_allowed,
          steal_index, task_color);
      }
      thread_depth[stealer] = saved_depth;
      if (stealer != tid) {
        #pragma omp atomic update
        idle_threads++;
      }
      g_available[tid][latest_level][idx[i]] = -1;
    }
  }
  #pragma omp taskwait
  return upper_bound;
}

int colorSort(vector<int> &candidate, vector<int> &orderOut,vector<int> &colorOut,
  const int tid, int* r_ed, int& stamp, const int end, const int start, int& pivot_vertex, vector<Bitset> &adj) {
  int total_cols = 0;
  stamp = 0;
  thread_local vector<vector<int>> color_class(1);
  pivot_vertex = candidate[end - 1];
  color_class.clear();
  color_class.resize(1);
  g_color_temp[tid][1] = 0;
  for (int i = end - 1; i >= start; i--) {
    const int v = candidate[i];
    r_ed[v] = 0;
    stamp++;

    for (int j = end - 1; j > i; j--) {
      const int u = candidate[j];
      if (adj[v].test(u)) {
        r_ed[v]++;
        r_ed[u]++;
        g_color_temp[tid][colorOut[u]] = stamp;
      }
    }

    int c = 1;
    while (g_color_temp[tid][c] == stamp) {
      c++;
    }
    colorOut[v] = c;
    g_color_temp[tid][colorOut[v] + 1] = 0;
    if (c > total_cols) {
      total_cols++;
      color_class.push_back(vector<int>());
    }
    color_class[c - 1].push_back(v);
  }

  for (int i = (int)color_class.size() - 1; i >= 0; i--) {
    for (int j = 0; j < (int)color_class[i].size(); j++) {
      orderOut.push_back(color_class[i][j]);
    }
  }
  return total_cols;
}

void expand(const int n, vector<int> &partial, vector<int> candidate, int latest_level, int thread, Bitset &candidBit, vector<Bitset> &adj) {
  #pragma omp atomic compare
  if (latest_level > max_size) {
    max_size = latest_level;
  }
  if (candidate.empty() || (int) candidate.size() + latest_level <= max_size) {
    return;
  }
  int stamp = 0;
  vector<int> orderLocal;
  int pivot_vertex;
  int maxColor = colorSort(candidate, orderLocal,
    color[thread][latest_level], thread, g_r_ed[thread][latest_level],
    stamp, (int)candidate.size(), 0, pivot_vertex, adj);
  if (maxColor + latest_level <= max_size) return;
  if (latest_level == 0) {
    vector<int> pos(n);
    for (int k = 0; k < (int)orderLocal.size(); k++) {
      pos[orderLocal[k]] = k;
    }
    #pragma omp parallel default(none) firstprivate(stamp) shared(idle_threads, pos, candidBit, prune, candidate, maxColor, g_color_temp, n, orderLocal, color, adj, max_size, partial, latest_level, thread, g_r_ed)
    {
      const int tid = omp_get_thread_num();
      #pragma omp for schedule(dynamic, 1)
      for (int i = 0; i < (int)orderLocal.size() - max_size; i++) {
        int v = orderLocal[i];
        stamp++;

        if (color[thread][latest_level][v] <= max_size) {
          continue;
        }

        vector<int> P = {v};
        Bitset v_common = adj[v];
        v_common &= candidBit;

        U64* data = v_common.a;
        for (int b = 0; b < W; b++) {
          U64 x = data[b];
          while (x) {
            int bit = __builtin_ctzll(x);
            int u = (b << 6) + bit;
            if (pos[u] <= i) {
              data[b] &= ~(1ULL << bit);
            }
            x &= (x - 1);
          }
        }

        vector<int> new_candidate;
        int color_ub = get_indices(v_common, new_candidate, tid, color[thread][latest_level], stamp);
        if (color_ub + latest_level + 1 > max_size)
        expand(n, P, move(new_candidate), latest_level + 1, tid, v_common, adj);
      }
    }
  } else {
    int deg = g_r_ed[thread][latest_level][pivot_vertex];
    if (deg == (int)candidate.size() - 1) {
      partial.push_back(pivot_vertex);
      candidBit.reset(pivot_vertex);
      candidate.pop_back();
      expand(n, partial, move(candidate), latest_level + 1, thread, candidBit, adj);
      partial.pop_back();
      return;
    }
    int size = (int)orderLocal.size();
    for (int i = 0; i < (int)orderLocal.size(); i++) {
      int v = orderLocal[i];
      candidBit.reset(v);
      if (latest_level + color[thread][latest_level][v] <= max_size) {
        break;
      }
      partial.push_back(v);
      stamp++;
      Bitset v_common = adj[v];
      v_common &= candidBit;
      vector<int> new_candidate;
      int color_ub = get_indices(v_common, new_candidate, thread, color[thread][latest_level], stamp);
      if (color_ub + latest_level + 1 > max_size)
        expand(n, partial, move(new_candidate), latest_level + 1, thread, v_common, adj);
      partial.pop_back();
    }
  }
}

int BRB(const ECLgraph g, int* const __restrict__ ed, int start_size, int& current_max, vector<int> &candidate,
  int valid_length, int* color_temp, bool task_allowed, vector<int> &colors,
  int latest_level, vector<int> &partial, const int tid) {
  for (int i = 0; i < start_size; i++) {
    if (prune[partial[i]]) return -1;
  }
  #pragma omp atomic compare
  if (start_size > max_size) {
    max_size = start_size;
  }
  current_max = max_size;
  if (valid_length == 0 || valid_length + start_size <= max_size) {
    return start_size;
  }

  int total_cols = 0;
  int color_array_size = std::min(valid_length, adj_size);
  memset(color_temp, 0, sizeof(int) * (color_array_size + 1));
  int stamp;
  int pivot = -1;
  int pivot_deg = -1;
  bool change = false;
  int induced_edges;
  do {
    stamp = 1;
    induced_edges = 0;
    for (int i = valid_length - 1; i >= 0; i--) {
      const int v = candidate[i];
      if (v == -1 || prune[v]) continue;
      const int v_beg = g.nindex[v];
      const int v_end = g.nindex[v + 1];
      int candid_pointer = valid_length - 1;
      for (int j = v_end - 1; j >= v_beg; j--) {
        const int u = g.nlist[j];
        if (u <= v) break;
        if (prune[u]) continue;
        while ((candid_pointer > i + 1) && (candidate[candid_pointer] > u || candidate[candid_pointer] == -1)) candid_pointer--;
        if (candidate[candid_pointer] == u) {
          induced_edges++;
          color_temp[colors[candid_pointer]] = stamp;
        }
      }
      int c = 1;
      while (color_temp[c] == stamp) ++c;
      colors[i] = c;
      total_cols = std::max(total_cols, c);
      stamp++;
    }
    
    stamp = -1;
    change = false;
    pivot = -1;
    pivot_deg = -1;
    pivot_deg = -1;
    change = false;
    for (int i = valid_length - 1; i >= 0; i--) {
      const int v = candidate[i];
      if (v == -1 || prune[v]) continue;
      const int v_beg = g.nindex[v];
      const int v_end = g.nindex[v + 1];
      int candid_pointer = valid_length - 1;
      int unique = 0;
      int deg = 0;
      for (int j = v_end - 1; j >= v_beg; j--) {
        const int u = g.nlist[j];
        if (prune[u]) continue;
        while ((candid_pointer > 0) && (candidate[candid_pointer] > u || candidate[candid_pointer] == -1)) candid_pointer--;
        if (candidate[candid_pointer] == u) {
          deg++;
          if (color_temp[colors[candid_pointer]] != stamp) {
            unique++;
          }
          color_temp[colors[candid_pointer]] = stamp;
        }
      }
      if (unique + start_size < max_size) {
        candidate[i] = -1;
        change = true;
      } else {
        if (deg > pivot_deg) {
          pivot_deg = deg;
          pivot = v;
        }
      }
      stamp--;
    }
  } while (change);
  if (total_cols + start_size < max_size || pivot == -1) return start_size;
  induced_edges *= 2;
  double density = (double) induced_edges / (valid_length * (valid_length - 1));
  //find cc
  int pivot_upper;
  vector<int> cc_color_size = component_count_csr(g, pivot,
    candidate, valid_length, colors, color_temp, prune, pivot_upper, color_array_size);


  int upper_bound = start_size;
  for (int i = -1; i < valid_length; i++) {
    int v;
    if (i != -1) {
      v = candidate[i];
      if (cc_color_size[i] + pivot_upper + 1 < max_size - start_size) {
        candidate[i] = -1;
        continue;
      }
      if ( v == -1 || prune[v] || (v == pivot)) continue;
    } else {
      v = pivot;
    }

    #pragma omp task if (idle_threads > 0) priority(latest_level) \
    shared(g, density, prune, max_size, candidate, valid_length, index_map, partial,\
    cc_color_size, colors, color_array_size, pivot_upper) \
    firstprivate(v, i, start_size, tid)
    {
      int stealer = omp_get_thread_num();
      int steal_index = stealer * adj_size;
      if (stealer != tid) {
        #pragma omp atomic update
	      idle_threads--;
      }
      bool leave = false;
      vector<int> new_partial;
      for (int j = 0; j < start_size; j++) {
	      new_partial.push_back(partial[j]);
        if (prune[partial[j]]) {leave = true; break;}
      }
      bool leave_task = false;
      if (!leave) {
        new_partial.push_back(v);
        memset(g_color_temp[stealer], 0, sizeof(int) * (color_array_size + 1));
        vector<int> next_candidate(valid_length);
        vector<int> next_color(valid_length);
        int task_stamp = -1;
        int next_length = 0;
        int next_ub = 1;
        const int v_beg = g.nindex[v];
        const int v_end = g.nindex[v + 1];
        int candid_pointer = 0;
        for (int j = v_beg; j < v_end; j++) {
          const int u = g.nlist[j];
          if (!prune[u]) {
            while ((candid_pointer < (valid_length - 1)) && (candidate[candid_pointer] < u || candidate[candid_pointer] == -1)) candid_pointer++;
            if (candidate[candid_pointer] == u) {
              if (u == pivot) {
                leave_task = true;
                break;
              }
              if (cc_color_size[candid_pointer] + pivot_upper + 1 <= max_size - start_size) {continue;}
              if (g_color_temp[stealer][colors[candid_pointer]] != task_stamp) {
                next_ub++;
              }
              g_color_temp[stealer][colors[candid_pointer]] = task_stamp;
              next_color[next_length] = colors[candid_pointer];
              next_candidate[next_length++] = u;
            }
          }
        }
        // call BRB on next_candidate with increased start size
        if (next_ub + start_size + 1 > max_size && !leave_task) {
          // check to see if without adjacency matrix I find the correct clique size
          // then check if without tasks
          if (next_length < adj_size && density > 0.5 && next_length > maxed && task_ready[stealer]) {
            task_ready[stealer] = false;
            int saved_depth = thread_depth[stealer];
            for (int j = 0; j < next_length; j++) {
              memset(word_Adj[steal_index + j], 0, sizeof(int) * next_length);
            }

            // fill out adjacency matrix
            for (int j = 0; j < next_length; j++) {
              const int u = next_candidate[j];
              index_map[steal_index + j] = u;
              g_available[stealer][saved_depth][j] = j + steal_index;
              word_Adj[steal_index + j][j] = 0;
              candid_pointer = 0;
              const int next_beg = g.nindex[u];
              const int next_end = g.nindex[u + 1];
              for (int k = next_beg; k < next_end; k++) {
                const int z = g.nlist[k];
                if (z > u) break;
                if (!prune[z]) {
                  while (candid_pointer < j - 1 && (next_candidate[candid_pointer] < z ||
                    next_candidate[candid_pointer] == -1)) candid_pointer++;
                  if (next_candidate[candid_pointer] == z) {
                    word_Adj[steal_index + j][candid_pointer] = 1;
                    word_Adj[steal_index + candid_pointer][j] = 1;
                  }
                }
              }
            }

            upper_bound = matrix_BRB(g, ed, start_size + 1, saved_depth, word_Adj, g_available[stealer],
              next_length, g_r_ed[stealer], max_size, steal_index, stealer, global_g_temp[stealer],
              g_color_temp[stealer],g_wl1[stealer], g_wl2[stealer],next_length >= max_wl, steal_index,
              next_color);
            task_ready[stealer] = true;
            thread_depth[stealer] = saved_depth;
          } else {
            upper_bound = BRB(g, ed, start_size + 1, max_size, next_candidate, next_length,
            color_temp, task_allowed, next_color, latest_level + 1, new_partial, stealer);
          }
        }
      }
      if (i != -1 && !leave_task) {
        candidate[i] = -1;
      }
      if (stealer != tid) {
        idle_threads++;
      }
    }
  }
  #pragma omp taskwait
  return upper_bound;
}

static int find(ECLgraph g, int* const __restrict__ ed, int* const __restrict__ sort,
  int* const __restrict__ count_array, int num_threads, int* overall_color, int** word_Adj, int* map) {
  global_g_temp = new ECLgraph[num_threads];
  for (int tid = 0; tid < num_threads; tid++) {
    global_g_temp[tid].eweight = NULL;
    global_g_temp[tid].nindex = (int*)malloc((adj_size + 1) * sizeof(int));
    global_g_temp[tid].nlist = (int*)malloc((adj_size * adj_size) * sizeof(int));
  }
  printf("num threads: %d\n", num_threads);

  //SHALLOW MSSB
  adj_size = MAXN;

  #pragma omp parallel default(none) shared(g, word_Adj, count_array, index_map, ed, max_size, adj_size)
  {
    int tid = omp_get_thread_num();
    int thread_index = adj_size * tid;
    #pragma omp for schedule(static, 128) reduction(max: max_size)
    for (int u = g.nodes - 1; u >= 0; u--) {
      if (ed[u] <= max_size) continue;
      const int u_beg = g.nindex[u];
      const int u_end = g.nindex[u + 1];
      const int u_deg = u_end - u_beg;
      int cur_length;
      int valid_length = std::min(adj_size, u_deg);
      memset(count_array + thread_index, 0, sizeof(int) * valid_length);
      for (int i = 0; i < valid_length; i++) {
        memset(word_Adj[thread_index + i], 0, sizeof(int) * valid_length);
      }
      int found_length = 0;
      for (int i = u_beg; i < u_end; i++) {
        const int v = g.nlist[i];
        if (found_length + 1 > adj_size) break;
        found_length++;
        cur_length = i - u_beg;
        const int v_beg = g.nindex[v];
        const int v_end = g.nindex[v + 1];
        const int v_deg = v_end - v_beg;
        index_map[thread_index + cur_length] = thread_index + cur_length;
        common_matrix(&g.nlist[u_beg], i - u_beg, &g.nlist[v_beg], v_deg, word_Adj, cur_length, count_array, thread_index);
        if (found_length == valid_length) break;
      }

      std::sort(index_map + thread_index, index_map + thread_index + valid_length, [&count_array](int a, int b) {
        return count_array[a] > count_array[b];
      });

      int bad_vertex = 0;
      bool bad = false;
      int clique_size = 1;
      for (int i = 1; i < valid_length; i++) {
        int neighbor_sum = 0;
        const int v_i = index_map[thread_index + i];
        for (int j = 0; j < i; j++) {
          const int v_j = index_map[thread_index + j] - thread_index;
          neighbor_sum += word_Adj[v_i][v_j];
        }
        if (neighbor_sum == i) {
          clique_size = i + 1;
        } else {
          if (!bad) {
            bad_vertex = index_map[thread_index + i];
          }
          if (bad && bad_vertex == index_map[thread_index + i]) {
            break;
          }
          int temp = index_map[thread_index + i];
          for (int k = i; k < valid_length - 1; k++) {
            index_map[thread_index + k] = index_map[thread_index + k + 1];
          }
          index_map[thread_index + valid_length - 1] = temp;
          i--;
          bad = true;
        }
      }
      #pragma omp atomic compare
      if (clique_size > max_size) {
        max_size = clique_size;
      }
      if (clique_size >= (valid_length - 1) || count_array[index_map[thread_index + clique_size + 1]] <= clique_size) {
        ed[u] = clique_size;
      }
    }
  }

  if (max_size + 1 == maxed) {
    #pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      freeECLgraph(global_g_temp[tid]);
    }
    delete [] global_g_temp;
    return max_size + 1;
  }

  vector<Global_task> reduced_g;
  int sort_order = 0;
  for (int v = 0; v < g.nodes; v++) {
    if (ed[v] <= max_size) {
      prune[v] = 1;
    }
    else {
      sort[v] = sort_order;
      sort_order++;
    }
  }

  printf("Number of nodes before pruning: %d\n", g.nodes);

  int edge_number = 0;
  int vertex_count = 0;
  int mdeg2 = 0;
  for (int v = 0; v < g.nodes; v++) {
    if (!prune[v]) {
      ed[vertex_count] = ed[v];
      overall_color[vertex_count] = overall_color[v];
      const int beg = g.nindex[v];
      const int end = g.nindex[v + 1];
      g.nindex[vertex_count] = edge_number;
      for (int i = beg; i < end; i++) {
        const int u = g.nlist[i];
        if (!prune[u]) {
          g.nlist[edge_number++] = sort[u];
        }
      }
      mdeg2 = std::max(mdeg2, edge_number - g.nindex[vertex_count]);
      reduced_g.push_back(Global_task(vertex_count, edge_number - g.nindex[vertex_count], ed[v]));
      vertex_count++;
    }
  }
  g.nindex[vertex_count] = edge_number;
  g.edges = edge_number;
  g.nodes = vertex_count;
  if (vertex_count == 0) return max_size + 1;
  memset(prune, 0, sizeof(int) * g.nodes);
  printf("After pruning: %d with max deg: %d\n", g.nodes, mdeg2);
  printf("Current estimate clique size: %d\n", max_size + 1);
  int reduced_index = reduced_g.size();
  int v_counter = 0;
  adj_size = MAXN;
  std::sort(reduced_g.begin(), reduced_g.begin() + reduced_index, [ed](const Global_task& a, const Global_task& b){
    return a.degree < b.degree;
  });

  #pragma omp parallel \
  shared(g_color_temp, idle_threads, max_wl, reduced_g, reduced_index, num_threads, \
  adj_size, v_counter, g, global_g_temp, count_array, index_map, ed, g_available, g_r_ed, max_size)
  {
    int tid = omp_get_thread_num();
    int thread_index = tid * adj_size;
    #pragma omp for schedule(dynamic, 1) nowait
    for (int c = 0; c < reduced_index - adj_size; c++) {
      int u = reduced_g[c].v;
      int check;
      #pragma omp atomic read
      check = max_size;
      thread_depth[tid] = 0;
      if (ed[u] > check) {
        bool task_allowed = false;
        const int u_beg = g.nindex[u];
        const int u_end = g.nindex[u + 1];
        const int u_deg = u_end - u_beg;
        if (u_deg >= max_wl) {
          task_allowed = true;
        }
        vector<int> color_vec(u_deg);
        vector<int> candidate_set(u_deg);
        int valid_length = 0;
        for (int i = u_beg; i < u_end; i++) {
          const int v = g.nlist[i];
          if (!prune[v]) {
            color_vec[valid_length] = overall_color[v];
            color_vec[valid_length] = 0;
            candidate_set[valid_length++] = v;
          }
        }
        vector<int> partial;
        int size = BRB(g, ed, 0, max_size, candidate_set, valid_length,
          g_color_temp[tid], task_allowed, color_vec, 0, partial, tid);
      }
      prune[u] = 1;
      int check_max;
      #pragma omp atomic read
      check_max = max_size;
      #pragma omp critical
      {
        v_counter++;
        std::cout << "\r\033[KPruned " << ((float) v_counter / reduced_index) * 100.0 << "% of vertices with size: " << check_max + 1 << std::flush;
      }
    }
    #pragma omp atomic update
    idle_threads++;
    #pragma omp barrier
  }

  sort_order = 0;
  printf("\nNumber of nodes before second pruning: %d\n", g.nodes);
  for (int v = 0; v < g.nodes; v++) {
    if (ed[v] <= max_size) {
      prune[v] = 1;
    }
    else if (!prune[v]) {
      sort[v] = sort_order;
      sort_order++;
    }
  }
  edge_number = 0;
  vertex_count = 0;
  mdeg2 = 0;
  for (int v = 0; v < g.nodes; v++) {
    if (!prune[v]) {
      ed[vertex_count] = ed[v];
      overall_color[vertex_count] = overall_color[v];
      const int beg = g.nindex[v];
      const int end = g.nindex[v + 1];
      g.nindex[vertex_count] = edge_number;
      for (int i = beg; i < end; i++) {
        const int u = g.nlist[i];
        if (!prune[u]) {
          g.nlist[edge_number++] = sort[u];
        }
      }
      mdeg2 = std::max(mdeg2, edge_number - g.nindex[vertex_count]);
      vertex_count++;
    }
  }
  g.nindex[vertex_count] = edge_number;
  g.edges = edge_number;
  g.nodes = vertex_count;
  memset(prune, 0, sizeof(int) * std::min(adj_size, g.nodes));
  printf("Number of nodes after second pruning: %d with max deg: %d\n", g.nodes, mdeg2);
  adj_size = MAXN;
  max_size += 1;
  vector<Bitset> adj(MAXN);
  vector<int> deg(g.nodes);
  for (int v = 0; v < g.nodes; v++) {
    map[v] = v;
    deg[v] = g.nindex[v + 1] - g.nindex[v];
  }
  vector<bool> removed(g.nodes, false);
  vector<int> core(g.nodes);
  int max_core = 0;
  int new_order = 0;
  while (new_order != g.nodes) {
    int min_index = -1;
    for (int v = 0; v < g.nodes; v++) {
      if (removed[v]) continue;
      if (min_index == -1 || deg[v] < deg[min_index]) min_index = v;
    }
    const int v_beg = g.nindex[min_index];
    const int v_end = g.nindex[min_index + 1];
    for (int j = v_beg; j < v_end; j++) {
      const int u = g.nlist[j];
      if (removed[u]) continue;
      deg[u]--;
    }
    core[min_index] = deg[min_index];
    if (core[min_index] > max_core) {max_core = core[min_index];}
    removed[min_index] = true;
    map[new_order++] = min_index;
  }
  vector<int> reOrderedCore(g.nodes);
  for (int i = 0; i < g.nodes; i++) {
    reOrderedCore[i] = core[map[i]];
  }
  for (int i = 0; i < g.nodes; i++) {
    sort[map[i]] = i;
  }

  ECLgraph new_g;
  new_g.nodes = g.nodes;
  new_g.nindex = new int [g.nodes + 1];
  new_g.nlist = new int [g.edges];
  new_g.edges = g.edges;
  new_g.eweight = NULL;
  vertex_count = 0;
  edge_number = 0;
  for (int v = 0; v < g.nodes; v++) {
    const int beg = g.nindex[map[v]];
    const int end = g.nindex[map[v] + 1];
    new_g.nindex[vertex_count] = edge_number;
    for (int i = beg; i < end; i++) {
      const int u = g.nlist[i];
      new_g.nlist[edge_number++] = sort[u];
    }
    vertex_count++;
    std::sort(new_g.nlist + new_g.nindex[v], new_g.nlist + edge_number);
  }
  new_g.nindex[vertex_count] = edge_number;

  vector<int> candidate;
  Bitset candidBit;
  for (int v = 0; v < g.nodes; v++) {
    candidate.push_back(v);
    candidBit.set(v);
    const int beg = new_g.nindex[v];
    const int end = new_g.nindex[v + 1];
    for (int i = beg; i < end; i++) {
      const int u = new_g.nlist[i];
      if (u > v) break;
      adj[v].set(u);
      adj[u].set(v);
    }
  }
  vector<int> partial;
  idle_threads = 0;
  expand(g.nodes, partial, move(candidate), 0, 0, candidBit, adj);

  #pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    freeECLgraph(global_g_temp[tid]);
  }
  printf("\n");
  freeECLgraph(new_g);
  delete [] global_g_temp;
  return max_size;
}

int main(int argc, char* argv []) {
  printf("CPU parallel maximum submatrix v1 (%s)\n", __FILE__);
  printf("Copyright 2026 Texas State University\n\n");

  if (argc != 2) {
    printf("USAGE: graph\n");
    exit(-1);
  }

  ECLgraph gr = readECLgraph(argv[1]);

  num_threads = 1;
  #pragma omp parallel shared(num_threads)
  {
    #pragma omp single
    {
      num_threads = omp_get_max_threads();
    }
  }

  //Allocation of memory
  std::vector<int> chunk_sum(num_threads, 0);
  std::vector<int> chunk_offset(num_threads, 0);
  std::vector<std::vector<int>> color_seen(omp_get_max_threads());
  for (auto& cs : color_seen)
    cs.assign(MAXN, 0);
  vector<int> color_iter(num_threads, 0);
  int* const sort = new int [gr.nodes];
  int* const map = new int [std::max(gr.nodes, gr.edges / 2)];
  int* const ed = new int [gr.nodes];
  int* g_deg = new int [gr.nodes];
  int* colors = new int [gr.nodes];
  int* const nlist2 = new int[gr.edges];
  int* const posscol = new int [gr.nodes];
  int* const posscol2 = new int [gr.edges / 2];
  int* const wl = new int [gr.nodes];

  adj_size = MAXN;
  g_available = new int** [num_threads];
  g_r_ed = new int** [num_threads];
  g_color_temp = new int* [num_threads];
  count_array = new int [adj_size * num_threads];
  index_map = new int [adj_size * num_threads];
  clique_embed = new int [MAXN + 1];
  g_wl1 = new int* [num_threads];
  g_wl2 = new int* [num_threads];
  task_ready = new bool [num_threads];
  prune = new int [gr.nodes];
  word_Adj = new int* [adj_size * num_threads];
  memset(prune, 0, sizeof(int) * gr.nodes);
  thread_depth = new int [num_threads];
  #pragma omp parallel for schedule(static, 128) shared(adj_size, num_threads)
  for (int i = 0; i < (adj_size * num_threads); i++) {
    index_map[i] = i;
    count_array[i] = 0;
    word_Adj[i] = new int [adj_size];
  }

  #pragma omp parallel for schedule(static, 1) shared(num_threads, adj_size)
  for (int tid = 0; tid < num_threads; tid++) {
    task_ready[tid] = true;
    g_color_temp[tid] = new int [adj_size + 1];
    g_wl1[tid] = new int [adj_size];
    g_wl2[tid] = new int [adj_size];
    g_available[tid] = new int* [MAXN];
    g_r_ed[tid] = new int* [MAXN];
    for (int step = 0; step < MAXN; step++) {
      g_available[tid][step] = new int [adj_size];
      g_r_ed[tid][step] = new int[adj_size];
    }
  }

  ECLgraph g;
  g.nodes = gr.nodes;
  g.edges = gr.edges;
  g.nindex = (int*)malloc((g.nodes + 1) * sizeof(int));
  g.nlist = (int*)malloc(g.edges * sizeof(int));
  g.eweight = NULL;
  g.nindex[0] = 0;
  printf("input: %s\n", argv[1]);
  printf("nodes: %d\n", gr.nodes);
  printf("edges: %d (%d)\n", gr.edges / 2, gr.edges);

  //touch every large array
  #pragma omp parallel for schedule(static)
  for (int u = 0; u < gr.nodes; u++) {
    sort[u] = 0;
    map[u] = 0;
    g.nindex[u] = 0;
    colors[u] = 0;
    g_deg[u] = gr.nindex[u + 1] - gr.nindex[u];
  }
  #pragma omp parallel for schedule(static)
  for (long e = 0; e < gr.edges; e++) {
    g.nlist[e] = 0;
    gr.nlist[e] = gr.nlist[e];
  }


  // start timed code sections
  timer.start();

  maxed = -1;
  colors = gc(gr, maxed, num_threads, colors, nlist2, posscol, posscol2, wl);

  #pragma omp parallel for schedule(static, 128)
  for (int u = 0; u < gr.nodes; u++) sort[u] = u;
  __gnu_parallel::sort(sort, sort + gr.nodes,
    [&](int a, int b)
    {return g_deg[a] > g_deg[b];},
    __gnu_parallel::multiway_mergesort_tag(num_threads));
  #pragma omp parallel for num_threads(num_threads) schedule(static)
  for (int u = 0; u < gr.nodes; u++) {
    const int q = sort[u];
    map[q] = u;
    g.nindex[u + 1] = gr.nindex[q + 1] - gr.nindex[q];
  }

  // make directed to only lower ID neighbors
  #pragma omp parallel num_threads(num_threads)
  {
    const int tid = omp_get_thread_num();
    const int chunk = (g.nodes + num_threads - 1) / num_threads;
    const int lo = std::min(tid * chunk, g.nodes);
    const int hi = std::min(lo + chunk, g.nodes);
    int local_sum = 0;
    for (int u = lo; u < hi; u++) {
      local_sum += g.nindex[u + 1];
      g.nindex[u + 1] = local_sum;
    }
    chunk_sum[tid] = local_sum;
    #pragma omp barrier
    #pragma omp single
    {
      int running = 0;
      for (int t = 0; t < num_threads; t++) {
        chunk_offset[t] = running;
        running += chunk_sum[t];
      }
    }
    for (int u = lo; u < hi; u++) {
      g.nindex[u + 1] += chunk_offset[tid];
    }
  }

  #pragma omp parallel for schedule(dynamic, 32)
  for (int u = 0; u < g.nodes; u++) {
    const int q = sort[u];
    const int q_beg = gr.nindex[q];
    const int q_end = gr.nindex[q + 1];
    int deg = q_end - q_beg;
    const int tid = omp_get_thread_num();
    int pos = g.nindex[u];
    int counter = 0;
    color_iter[tid]++;
    std::vector<int>& cseen = color_seen[tid];
    for (int i = q_beg; i < q_end; i++) {
      const int v = map[gr.nlist[i]];
      const int c = colors[gr.nlist[i]];
      counter += (cseen[c] != color_iter[tid]);
      cseen[c] = color_iter[tid];
      g.nlist[pos] = v;
      pos++;
    }
    ed[u] = counter;
    std::sort(&g.nlist[g.nindex[u]], &g.nlist[pos]);
  }
  
  float time_so_far = timer.elapsed();
  printf("gr+ed time: %.6f s\n", time_so_far);
  printf("max ed:  %5d\n", maxed);
  max_wl = maxed;
  max_size = 0;
  timer.start();
  max_size = find(g, ed, sort, count_array, num_threads, colors, word_Adj, map);
  float search_time = timer.elapsed();
  printf("max_find size: %d\n", max_size);
  printf("Search time: %.4f\n", search_time);
  printf("Total time: %.4f\n", search_time + time_so_far);

  for (int i = (adj_size * num_threads) - 1; i >= 0; i--) {
    delete[] word_Adj[i];
  }

  #pragma omp parallel for schedule(static, 1)
  for (int tid = 0; tid < num_threads; tid++) {
    delete [] g_wl1[tid];
    delete [] g_wl2[tid];
    delete [] g_color_temp[tid];
    for (int step = 0; step < MAXN; step++) {
      delete [] g_available[tid][step];
      delete [] g_r_ed[tid][step];
    }
    delete [] g_available[tid];
    delete [] g_r_ed[tid];
  }
  delete [] nlist2;
  delete [] posscol;
  delete [] posscol2;
  delete [] wl;
  delete [] thread_depth;
  delete [] task_ready;
  delete [] prune;
  delete [] g_color_temp;
  delete [] g_available;
  delete [] g_r_ed;
  delete [] g_wl1;
  delete [] g_wl2;
  delete [] colors;
  delete [] word_Adj;
  delete [] count_array;
  delete [] clique_embed;
  delete [] index_map;
  delete [] sort;
  delete [] map;
  freeECLgraph(g);
  delete [] ed;
  delete [] g_deg;
  return 0;
}











