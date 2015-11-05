#include <malloc.h>
#include <stdint.h>
#include <xmmintrin.h>
#include <sys/time.h>
#include <climits>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <queue>
#include <set>
#include <algorithm>
#include <fstream>
#include <utility>
#include <limits>
#include <omp.h>
#include "agl.h"

// Return the number of threads that would be executed in parallel regions
int get_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

// Set the number of threads that would be executed in parallel regions
void set_num_threads(int num_threads) {
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#else
  if (num_threads != 1) {
    assert(!"compile with -fopenmp");
  }
#endif
}

int get_thread_id() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

template <typename T>
struct parallel_vector {
  parallel_vector(int max_threads, size_t size_limit) : v(max_threads, std::vector<T>(size_limit)), n(max_threads, 0) {}

  void push_back(const T &x) {
    int id = get_thread_id();
    v[id][n[id]++] = x;
  }

  void clear() {
    for (int i = 0; i < n.size(); ++i) n[i] = 0;
  }

  std::vector<std::vector<T> > v;
  std::vector<size_t> n;
};

template <int kNumBitParallelRoots = 50>
class PrunedLandmarkLabeling {
 public:
  // Constructs an index from a graph, given as a list of edges.
  // Vertices should be described by numbers starting from zero.
  // Returns |true| when successful.
  bool ConstructIndex(const std::vector<std::pair<int, int> > &es);
  bool ConstructIndex(std::istream &ifs);
  bool ConstructIndex(const char *filename);

  // Returns distance vetween vertices |v| and |w| if they are connected.
  // Otherwise, returns |INT_MAX|.
  inline int QueryDistance(int v, int w);

  // Insert an edge (v, w)
  bool InsertEdge(int v, int w);

  int GetNumVertices() { return num_v_; }
  void Free();
  double GetAverageLabelSize();

  PrunedLandmarkLabeling() : num_v_(0), index_(NULL), index_l_(0), time_load_(0), time_indexing_(0) {}
  virtual ~PrunedLandmarkLabeling() { Free(); }

  // private:
  static const uint8_t INF8;    // For unreachable pairs
  static const uint32_t INF32;  // For sentinel
  static const int kInitialLabelCapacity;

  struct index_t {
    uint8_t bpspt_d[kNumBitParallelRoots];
    uint64_t bpspt_s[kNumBitParallelRoots][2];  // [0]: S^{-1}, [1]: S^{0}
    uint32_t *spt_v;
    uint8_t *spt_d;
    uint32_t spt_l;

    void Expand() {
      int new_spt_l = std::max(kInitialLabelCapacity, int((spt_l + 1) * 1.5));
      uint32_t *new_spt_v = (uint32_t *)memalign(64, new_spt_l * sizeof(uint32_t));
      uint8_t *new_spt_d = (uint8_t *)memalign(64, new_spt_l * sizeof(uint8_t));
      assert(new_spt_v && new_spt_d);
      memcpy(new_spt_v, spt_v, spt_l * sizeof(int32_t));
      memcpy(new_spt_d, spt_d, spt_l * sizeof(int8_t));
      memset(new_spt_v + spt_l, 0, (new_spt_l - spt_l) * sizeof(int32_t));
      free(spt_v);
      free(spt_d);
      spt_v = new_spt_v;
      spt_d = new_spt_d;
      // printf(" EXPAND: %d -> %d\n", spt_l, new_spt_l);
      spt_l = new_spt_l;
    }
  } __attribute__((aligned(64)));  // Aligned for cache lines

  int num_v_;
  std::vector<std::vector<int> > adj_;
  index_t *index_;
  int index_l_;           // We don't need to use |size_t|
  std::vector<int> ord_;  // ord_[i] = v (e.g., ord[spt_v[0]] = the vertex with the highest degree)

  void PartialBFS(uint32_t r, int sv, int sd);
  void PartialBPBFS(int k, int v, int w);

  // Statistics
  double time_load_, time_indexing_;

  double GetCurrentTimeSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
  }
};

template <int kNumBitParallelRoots>
const uint8_t PrunedLandmarkLabeling<kNumBitParallelRoots>::INF8 = 100;  // We require |INF8 + INF8| fits in |uint8_t|
template <int kNumBitParallelRoots>
const uint32_t PrunedLandmarkLabeling<kNumBitParallelRoots>::INF32 =
    std::numeric_limits<int32_t>::max();  // signed for safety
template <int kNumBitParallelRoots>
const int PrunedLandmarkLabeling<kNumBitParallelRoots>::kInitialLabelCapacity = 8;

template <int kNumBitParallelRoots>
bool PrunedLandmarkLabeling<kNumBitParallelRoots>::ConstructIndex(const char *filename) {
  std::ifstream ifs(filename);
  return ifs && ConstructIndex(ifs);
}

template <int kNumBitParallelRoots>
bool PrunedLandmarkLabeling<kNumBitParallelRoots>::ConstructIndex(std::istream &ifs) {
  std::vector<std::pair<int, int> > es;
  for (int v, w; ifs >> v >> w;) {
    es.push_back(std::make_pair(v, w));
  }
  if (ifs.bad()) return false;
  ConstructIndex(es);
  return true;
}

template <int kNumBitParallelRoots>
bool PrunedLandmarkLabeling<kNumBitParallelRoots>::ConstructIndex(const std::vector<std::pair<int, int> > &es) {
  //
  // Prepare the adjacency list and index space
  //
  Free();
  time_load_ = -GetCurrentTimeSec();
  int E = es.size();
  int &V = num_v_;
  V = 0;
  for (size_t i = 0; i < es.size(); ++i) {
    V = std::max(V, std::max(es[i].first, es[i].second) + 1);
  }
  adj_.resize(V);
  for (size_t i = 0; i < es.size(); ++i) {
    int v = es[i].first, w = es[i].second;
    adj_[v].push_back(w);
    adj_[w].push_back(v);
  }
  time_load_ += GetCurrentTimeSec();

  index_ = (index_t *)memalign(64, V * sizeof(index_t));
  index_l_ = V;
  if (index_ == NULL) {
    num_v_ = 0;
    return false;
  }
  for (int v = 0; v < V; ++v) {
    index_[v].spt_v = NULL;
    index_[v].spt_d = NULL;
    index_[v].spt_l = 0;
  }

  //
  // Order vertices by decreasing order of degree
  //
  time_indexing_ = -GetCurrentTimeSec();
  std::vector<int> &inv = ord_;           // new label -> old label
  std::vector<std::vector<int> > adj(V);  // Do not confuse with |adj_| (TODO: fix)
  {
    // Order
    std::vector<std::pair<float, int> > deg(V);
    for (int v = 0; v < V; ++v) {
      // We add a random value here to diffuse nearby vertices
      deg[v] = std::make_pair(adj_[v].size() + float(rand()) / RAND_MAX, v);
    }
    // puts("TODO: SORT");
    std::sort(deg.rbegin(), deg.rend());
    for (int i = 0; i < V; ++i) {
      int v = deg[i].second;
      if (!adj_[v].empty()) inv.push_back(v);  // We ignore isolated vertices here (to decide the ordering later)
    }

    // Relabel the vertex IDs
    std::vector<int> rank(V);
    for (int i = 0; i < V; ++i) rank[deg[i].second] = i;
    for (int v = 0; v < V; ++v) {
      for (size_t i = 0; i < adj_[v].size(); ++i) {
        adj[rank[v]].push_back(rank[adj_[v][i]]);
      }
    }
  }

  //
  // Bit-parallel labeling
  //
  std::vector<bool> usd(V, false);  // Used as root? (in new label)
  {
    std::vector<uint8_t> tmp_d(V);
    std::vector<std::pair<uint64_t, uint64_t> > tmp_s(V);
    std::vector<int> que(V);
    std::vector<std::pair<int, int> > sibling_es(E);
    std::vector<std::pair<int, int> > child_es(E);

    int r = 0;
    for (int i_bpspt = 0; i_bpspt < kNumBitParallelRoots; ++i_bpspt) {
      while (r < V && usd[r]) ++r;
      if (r == V) {
        for (int v = 0; v < V; ++v) index_[v].bpspt_d[i_bpspt] = INF8;
        continue;
      }
      usd[r] = true;

      fill(tmp_d.begin(), tmp_d.end(), INF8);
      fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));

      int que_t0 = 0, que_t1 = 0, que_h = 0;
      que[que_h++] = r;
      tmp_d[r] = 0;
      que_t1 = que_h;

      int ns = 0;
      std::vector<int> vs;
      sort(adj[r].begin(), adj[r].end());
      for (size_t i = 0; i < adj[r].size(); ++i) {
        int v = adj[r][i];
        if (!usd[v]) {
          usd[v] = true;
          que[que_h++] = v;
          tmp_d[v] = 1;
          tmp_s[v].first = 1ULL << ns;
          vs.push_back(v);
          if (++ns == 64) break;
        }
      }

      for (int d = 0; que_t0 < que_h; ++d) {
        int num_sibling_es = 0, num_child_es = 0;

        for (int que_i = que_t0; que_i < que_t1; ++que_i) {
          int v = que[que_i];

          for (size_t i = 0; i < adj[v].size(); ++i) {
            int tv = adj[v][i];
            int td = d + 1;

            if (d > tmp_d[tv])
              ;
            else if (d == tmp_d[tv]) {
              if (v < tv) {
                sibling_es[num_sibling_es].first = v;
                sibling_es[num_sibling_es].second = tv;
                ++num_sibling_es;
              }
            } else {
              if (tmp_d[tv] == INF8) {
                que[que_h++] = tv;
                tmp_d[tv] = td;
              }
              child_es[num_child_es].first = v;
              child_es[num_child_es].second = tv;
              ++num_child_es;
            }
          }
        }

        for (int i = 0; i < num_sibling_es; ++i) {
          int v = sibling_es[i].first, w = sibling_es[i].second;
          tmp_s[v].second |= tmp_s[w].first;
          tmp_s[w].second |= tmp_s[v].first;
        }
        for (int i = 0; i < num_child_es; ++i) {
          int v = child_es[i].first, c = child_es[i].second;
          tmp_s[c].first |= tmp_s[v].first;
          tmp_s[c].second |= tmp_s[v].second;
        }

        que_t0 = que_t1;
        que_t1 = que_h;
      }

      for (int v = 0; v < (int)inv.size(); ++v) {
        index_[inv[v]].bpspt_d[i_bpspt] = tmp_d[v];
        index_[inv[v]].bpspt_s[i_bpspt][0] = tmp_s[v].first;
        index_[inv[v]].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
      }
      for (int v = 0; v < V; ++v)
        if (adj_[v].empty()) {
          index_[v].bpspt_d[i_bpspt] = INF8;
          index_[v].bpspt_s[i_bpspt][0] = index_[v].bpspt_s[i_bpspt][1] = 0;
        }
    }
  }

  //
  // Pruned labeling
  //
  int max_threads = get_max_threads();
  {
    // Sentinel (V, INF8) is added to all the vertices
    std::vector<std::pair<std::vector<int>, std::vector<uint8_t> > > tmp_idx(
        V, make_pair(std::vector<int>(1, INF32), std::vector<uint8_t>(1, INF8)));

    // std::vector<bool> vis(V);
    std::vector<uint8_t> vis(V / 8 + 1);
    std::vector<int> que(V);
    std::vector<uint8_t> dst_r(V + 1, INF8);

    parallel_vector<int> pdiff_nxt_que(max_threads, V);

    for (int r = 0; r < V; ++r) {
      if (usd[r] || adj[r].empty()) continue;
      index_t &idx_r = index_[inv[r]];
      const std::pair<std::vector<int>, std::vector<uint8_t> > &tmp_idx_r = tmp_idx[r];
      for (size_t i = 0; i < tmp_idx_r.first.size() - 1; ++i) {
        dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
      }

      int que_t0 = 0, que_t1 = 0, que_h = 0;
      que[que_h++] = r;
      // vis[r] = true;
      vis[r / 8] |= 1 << (r % 8);
      que_t1 = que_h;

      for (int d = 0; que_t0 < que_h; ++d) {
        for (int que_i = que_t0; que_i < que_t1; ++que_i) {
          int v = que[que_i];
          std::pair<std::vector<int>, std::vector<uint8_t> > &tmp_idx_v = tmp_idx[v];
          assert(v < (int)inv.size());
          index_t &idx_v = index_[inv[v]];

          // Prefetch
          _mm_prefetch(&idx_v.bpspt_d[0], _MM_HINT_T0);
          _mm_prefetch(&idx_v.bpspt_s[0][0], _MM_HINT_T0);
          _mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
          _mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);

          // Prune?
          if (usd[v]) continue;
          for (int i = 0; i < kNumBitParallelRoots; ++i) {
            int td = idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
            if (td - 2 <= d) {
              td += (idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 : ((idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
                                                                        (idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
                                                                           ? -1
                                                                           : 0;
              if (td <= d) goto pruned;
            }
          }
          for (size_t i = 0; i < tmp_idx_v.first.size() - 1; ++i) {
            int w = tmp_idx_v.first[i];
            int td = tmp_idx_v.second[i] + dst_r[w];
            if (td <= d) goto pruned;
          }

          // Traverse
          tmp_idx_v.first.back() = r;
          tmp_idx_v.second.back() = d;
          tmp_idx_v.first.push_back(INF32);
          tmp_idx_v.second.push_back(INF8);

          for (size_t i = 0; i < adj[v].size(); ++i) {
            int w = adj[v][i];
            // if (!vis[w]) {
            //   que[que_h++] = w;
            //   vis[w] = true;
            // }
            for (;;) {
              uint8_t m = vis[w / 8];
              if (m & (1 << (w % 8))) break;
              if (__sync_bool_compare_and_swap(&vis[w / 8], m, m | (1 << (w % 8)))) {
                pdiff_nxt_que.push_back(w);
                break;
              }
            }
          }
        pruned : {}
        }

        for (int i = 0; i < max_threads; ++i) {
          for (int j = 0; j < (int)pdiff_nxt_que.n[i]; ++j) {
            que[que_h++] = pdiff_nxt_que.v[i][j];
          }
        }

        que_t0 = que_t1;
        que_t1 = que_h;
        pdiff_nxt_que.clear();
      }

      for (int i = 0; i < que_h; ++i) vis[que[i] / 8] = 0;
      for (size_t i = 0; i < tmp_idx_r.first.size() - 1; ++i) {
        dst_r[tmp_idx_r.first[i]] = INF8;
      }
      usd[r] = true;
    }

    for (int v = 0; v < (int)inv.size(); ++v) {
      int k = tmp_idx[v].first.size();
      index_[inv[v]].spt_v = (uint32_t *)memalign(64, k * sizeof(uint32_t));
      index_[inv[v]].spt_d = (uint8_t *)memalign(64, k * sizeof(uint8_t));
      index_[inv[v]].spt_l = k;
      if (!index_[inv[v]].spt_v || !index_[inv[v]].spt_d) {
        Free();
        return false;
      }
      for (int i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
      for (int i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
      tmp_idx[v].first.clear();
      tmp_idx[v].second.clear();
    }
    for (int i = 0; i < V; ++i)
      if (adj_[i].empty()) {
        assert(index_[i].spt_v == NULL);
        index_[i].spt_v = (uint32_t *)memalign(64, 1 * sizeof(uint32_t));
        index_[i].spt_d = (uint8_t *)memalign(64, 1 * sizeof(uint8_t));
        index_[i].spt_v[0] = INF32;
      }
  }

  time_indexing_ += GetCurrentTimeSec();
  return true;
}

template <int kNumBitParallelRoots>
int PrunedLandmarkLabeling<kNumBitParallelRoots>::QueryDistance(int v, int w) {
  if (v == w) return 0;
  if (v >= num_v_ || w >= num_v_) return INT_MAX;

  const index_t &idx_v = index_[v];
  const index_t &idx_w = index_[w];
  int d = INF8;

  _mm_prefetch(&idx_v.spt_v[0], _MM_HINT_T0);
  _mm_prefetch(&idx_w.spt_v[0], _MM_HINT_T0);
  _mm_prefetch(&idx_v.spt_d[0], _MM_HINT_T0);
  _mm_prefetch(&idx_w.spt_d[0], _MM_HINT_T0);

  for (int i = 0; i < kNumBitParallelRoots; ++i) {
    int td = idx_v.bpspt_d[i] + idx_w.bpspt_d[i];
    if (td - 2 <= d) {
      td += (idx_v.bpspt_s[i][0] & idx_w.bpspt_s[i][0])
                ? -2
                : ((idx_v.bpspt_s[i][0] & idx_w.bpspt_s[i][1]) | (idx_v.bpspt_s[i][1] & idx_w.bpspt_s[i][0])) ? -1 : 0;

      if (td < d) d = td;
    }
  }
  for (int i1 = 0, i2 = 0;;) {
    uint32_t v1 = idx_v.spt_v[i1], v2 = idx_w.spt_v[i2];
    if (v1 == v2) {
      if (v1 == INF32) break;  // Sentinel
      // printf(" %d->%d->%d: %d+%d\n", v, ord_[v1], w, idx_v.spt_d[i1], idx_w.spt_d[i2]);
      int td = idx_v.spt_d[i1] + idx_w.spt_d[i2];
      if (td < d) d = td;
      ++i1;
      ++i2;
    } else {
      i1 += v1 < v2 ? 1 : 0;
      i2 += v1 > v2 ? 1 : 0;
    }
  }

  if (d >= INF8 - 2) d = INT_MAX;
  return d;
}

template <int kNumBitParallelRoots>
void PrunedLandmarkLabeling<kNumBitParallelRoots>::Free() {
  for (int v = 0; v < num_v_; ++v) {
    free(index_[v].spt_v);
    free(index_[v].spt_d);
  }
  free(index_);
  index_ = NULL;
  num_v_ = 0;
}

template <int kNumBitParallelRoots>
double PrunedLandmarkLabeling<kNumBitParallelRoots>::GetAverageLabelSize() {
  // std::cout << "load time: "     << time_load_     << " seconds" << std::endl;
  // std::cout << "indexing time: " << time_indexing_ << " seconds" << std::endl;

  double s = 0.0;
  for (int v = 0; v < num_v_; ++v) {
    for (int i = 0; index_[v].spt_v[i] != INF32; ++i) {
      ++s;
    }
  }
  return s /= num_v_;
}

template <int kNumBitParallelRoots>
void PrunedLandmarkLabeling<kNumBitParallelRoots>::PartialBPBFS(int r, int u, int v) {

  static std::vector<int> que0, que1;
  static std::vector<bool> queued;
  if ((int)queued.size() < num_v_) {
    que0.resize(num_v_ * 2);
    que1.resize(num_v_ * 2);
    queued.resize(num_v_ * 2);
  }

  int h0 = 0, h1 = 0;
  int base_d = std::min(index_[u].bpspt_d[r], index_[v].bpspt_d[r]);
  if (base_d == INF8) return;
  if (index_[u].bpspt_d[r] == base_d) queued[que0[h0++] = u] = true;
  if (index_[v].bpspt_d[r] == base_d) queued[que0[h0++] = v] = true;

  for (int d = base_d; h0 > 0; ++d) {
    for (int i0 = 0; i0 < h0; ++i0) {
      int v = que0[i0];
      // printf(" %d: %d\n", d, v);

      for (int i = 0; i < adj_[v].size(); ++i) {
        int tv = adj_[v][i], td = d + 1;
        if (d == index_[tv].bpspt_d[r]) {
          // Set propagation (1)
          uint64_t ts = index_[tv].bpspt_s[r][1] | index_[v].bpspt_s[r][0];
          if (ts != index_[tv].bpspt_s[r][1]) {
            index_[tv].bpspt_s[r][1] = ts;
            // if (!queued[tv]) queued[que0[h0++] = tv] = true;  // Update -> Enque
          }
        }
        if (td < index_[tv].bpspt_d[r]) {
          index_[tv].bpspt_s[r][1] = (index_[tv].bpspt_d[r] == td + 1 ? index_[tv].bpspt_s[r][0] : 0);
          index_[tv].bpspt_s[r][0] = 0;
          index_[tv].bpspt_d[r] = td;
          assert(!queued[tv]);
          queued[que1[h1++] = tv] = true;
        }
      }
    }

    for (int i0 = 0; i0 < h0; ++i0) {
      int v = que0[i0];
      for (int i = 0; i < adj_[v].size(); ++i) {
        int tv = adj_[v][i];
        if (index_[tv].bpspt_d[r] == d + 1) {
          // Set propagation (2)
          uint64_t ts0 = index_[tv].bpspt_s[r][0] | index_[v].bpspt_s[r][0];
          uint64_t ts1 = index_[tv].bpspt_s[r][1] | index_[v].bpspt_s[r][1];
          if (ts0 != index_[tv].bpspt_s[r][0] || ts1 != index_[tv].bpspt_s[r][1]) {
            index_[tv].bpspt_s[r][0] = ts0;
            index_[tv].bpspt_s[r][1] = ts1;
            // if (!queued[tv]) queued[que1[h1++] = tv] = true;  // Update -> Enque
          }
        }
      }
    }

    for (int i0 = 0; i0 < h0; ++i0) queued[que0[i0]] = false;
    que0.swap(que1);
    h0 = h1;
    que1.clear();
    h1 = 0;
  }
}

template <int kNumBitParallelRoots>
void PrunedLandmarkLabeling<kNumBitParallelRoots>::PartialBFS(uint32_t bfs_i, int sv, int sd) {

  static std::vector<std::pair<int, int> > que;
  static std::vector<bool> vis;
  static std::vector<int> root_label;
  if ((int)que.size() < num_v_) {
    que.resize(num_v_ * 2);
    vis.resize(num_v_ * 2);
    root_label.resize(num_v_ * 2);
  }

  int r = ord_[bfs_i];
  index_t &idx_r = index_[r];
  for (int i = 0; idx_r.spt_v[i] != INF32; ++i) {
    root_label[idx_r.spt_v[i]] = idx_r.spt_d[i];
  }

  int que_h = 0, que_t = 0;
  // queue<pair<int, int> > que;
  que[que_t++] = std::make_pair(sv, sd);
  vis[sv] = true;

  while (que_h < que_t) {
    int v = que[que_h].first;
    int d = que[que_h].second;
    ++que_h;

    // Puruning test & new label
    {
      index_t &idx_v = index_[v];

      for (int i = 0; i < kNumBitParallelRoots; ++i) {
        int td = idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
        if (td - 2 <= d) {
          td += (idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0])
                    ? -2
                    : ((idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) | (idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0])) ? -1
                                                                                                                  : 0;
          if (td <= d) goto prune;
        }
      }

      // Case 1: |bfs_i| is already in |label[v]|
      int i = 0;
      for (; idx_v.spt_v[i] <= bfs_i; ++i) {
        uint32_t li = idx_v.spt_v[i];
        int ld = idx_v.spt_d[i];

        if (li == bfs_i) {
          if (ld <= d)
            goto prune;
          else {
            idx_v.spt_d[i] = d;
            goto traverse;
          }
        }
        if (root_label[li] != -1 && root_label[li] + ld <= d) goto prune;
      }

      // Case 2: |bfs_i| is not present in |label[v]|
      if (idx_v.spt_v[idx_v.spt_l - 1] == INF32) idx_v.Expand();
      int j;
      for (j = i; idx_v.spt_v[j] != INF32; ++j)
        ;
      for (; j >= i; --j) {
        idx_v.spt_v[j + 1] = idx_v.spt_v[j];
        idx_v.spt_d[j + 1] = idx_v.spt_d[j];
      }
      idx_v.spt_v[i] = bfs_i;
      idx_v.spt_d[i] = d;
    }

  traverse:
    ;
    for (int i = 0; i < adj_[v].size(); ++i) {
      int w = adj_[v][i];
      if (!vis[w]) {
        que[que_t++] = std::make_pair(w, d + 1);
        vis[w] = true;
      }
    }

  prune:
    ;
  }

  for (int i = 0; i < que_t; ++i) vis[que[i].first] = false;
  for (int i = 0; index_[r].spt_v[i] != INF32; ++i) {
    root_label[index_[r].spt_v[i]] = -1;
  }
}

template <int kNumBitParallelRoots>
bool PrunedLandmarkLabeling<kNumBitParallelRoots>::InsertEdge(int u, int v) {
  int new_num_v = std::max(num_v_, std::max(u, v) + 1);
  if (new_num_v > index_l_) {
    int new_index_l = std::max(new_num_v, index_l_ * 2);
    index_t *new_index = (index_t *)memalign(64, new_index_l * sizeof(index_t));
    if (new_index == NULL) return false;
    memcpy(new_index, index_, index_l_ * sizeof(index_t));
    free(index_);
    index_ = new_index;
    // printf("EXPANDER: %d -> %d\n", index_l_, new_index_l);
    index_l_ = new_index_l;
  }
  for (int i = num_v_; i < new_num_v; ++i) {
    index_t &idx = index_[i];
    for (int k = 0; k < kNumBitParallelRoots; ++k) {
      idx.bpspt_d[k] = INF8;
      idx.bpspt_s[k][0] = idx.bpspt_s[k][1] = 0;
    }
    idx.spt_v = (uint32_t *)memalign(64, kInitialLabelCapacity * sizeof(uint32_t));
    idx.spt_d = (uint8_t *)memalign(64, kInitialLabelCapacity * sizeof(uint8_t));
    if (!idx.spt_v || !idx.spt_d) return false;
    idx.spt_l = kInitialLabelCapacity;
    memset(idx.spt_v, 0, kInitialLabelCapacity * sizeof(uint32_t));
    idx.spt_v[0] = INF32;
  }
  num_v_ = new_num_v;
  adj_.resize(num_v_);
  adj_[u].push_back(v);
  adj_[v].push_back(u);

  for (int i = 0; i < 2; ++i) {
    int x = i == 0 ? u : v;
    index_t &idx = index_[x];
    if (adj_[x].size() == 1) {  // New comer
      // printf("New comer: %d\n", x);
      int id = ord_.size();
      ord_.push_back(x);
      if (idx.spt_l <= 1) idx.Expand();
      idx.spt_v[0] = id;
      idx.spt_d[0] = 0;
      idx.spt_v[1] = INF32;
    }
  }

  // Update bit-parallel labels
  // puts("BIT PARALLEL");
  for (int k = 0; k < kNumBitParallelRoots; ++k) {
    PartialBPBFS(k, u, v);
  }

  // Update normal labels
  // puts("NORMAL");
  const index_t &idx_u = index_[u], &idx_v = index_[v];
  int iu = 0, iv = 0;
  for (;;) {
    uint32_t vu = idx_u.spt_v[iu], vv = idx_v.spt_v[iv];
    int du = idx_u.spt_d[iu], dv = idx_v.spt_d[iv];

    if (vu < vv) {  // u -> v
      PartialBFS(vu, v, du + 1);
      ++iu;
    } else if (vu > vv) {  // v -> u
      PartialBFS(vv, u, dv + 1);
      ++iv;
    } else {  // u <-> v
      if (vu == INF32) break;
      if (du + 1 < dv) PartialBFS(vu, v, du + 1);
      if (dv + 1 < du) PartialBFS(vv, u, dv + 1);
      ++iu;
      ++iv;
    }
  }

  return true;
}
