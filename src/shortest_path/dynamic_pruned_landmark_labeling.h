#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"
#include <malloc.h>

namespace agl {

template<typename T>
struct parallel_vector {
  parallel_vector(int max_threads, size_t size_limit)
      : v(max_threads, std::vector<T>(size_limit)), n(max_threads, 0) {}

  void push_back(const T &x) {
    // int id = get_thread_id();
    int id = 0;
    v[id][n[id]++] = x;
  }

  void clear() {
    for (int i = 0; i < n.size(); ++i) n[i] = 0;
  }

  std::vector<std::vector<T>> v;
  std::vector<size_t> n;
};

template<size_t kNumBitParallelRoots = 16>
class dynamic_pruned_landmark_labeling
    : public dynamic_graph_index_interface<G>,
      public distance_query_interface<G> {
 public:
  dynamic_pruned_landmark_labeling() : num_v_(0) {}
  virtual ~dynamic_pruned_landmark_labeling() { free_all(); }
  virtual void construct(const G &g) override;
  virtual W query_distance(const G &g, V v_from, V v_to) override;

  virtual void add_edge(const G &g, V v_from, const E &e) override;

  virtual void remove_edge(const G &g, V v_from, V v_to) override {
    assert(false);
  }

  virtual void add_vertices(const G &g, V old_num_vertices) override {
    assert(false);
  }

  virtual void remove_vertices(const G &, V old_num_vertices) override {
    assert(false);
  }

 private:
  static const uint8_t INF8 = 100;  // For unreachable pairs
  static const uint32_t INF32 =
      std::numeric_limits<int32_t>::max();  // For sentinel

  struct index_t {
    uint32_t *spt_v;
    uint8_t *spt_d;
    uint32_t spt_l;

    void Expand() {
      int new_spt_l = int((spt_l + 1) * 1.5);
      uint32_t *new_spt_v =
          (uint32_t *)memalign(64, new_spt_l * sizeof(uint32_t));
      uint8_t *new_spt_d = (uint8_t *)memalign(64, new_spt_l * sizeof(uint8_t));
      assert(new_spt_v && new_spt_d);
      memcpy(new_spt_v, spt_v, spt_l * sizeof(int32_t));
      memcpy(new_spt_d, spt_d, spt_l * sizeof(int8_t));
      memset(new_spt_v + spt_l, 0, (new_spt_l - spt_l) * sizeof(int32_t));
      free(spt_v);
      free(spt_d);
      spt_v = new_spt_v;
      spt_d = new_spt_d;
      spt_l = new_spt_l;
    }
  };

  V num_v_;
  std::vector<std::vector<V>> adj_[2];
  std::vector<index_t> idx_[2];
  std::vector<V> ord_;

  void free_all() {
    for (int i = 0; i < 2; ++i) {
      for (V v = 0; v < num_v_; ++v) {
        free(idx_[i][v].spt_v);
        free(idx_[i][v].spt_d);
      }
      idx_[i].clear();
    }
    num_v_ = 0;
  }

  void partial_bfs(V bfs_i, V sv, W sd, int x);
};

template<size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>
::construct(const G &g) {
  free_all();
  V &num_v = num_v_;
  num_v = g.num_vertices();
  const std::vector<std::pair<V, E>> &es = g.edge_list();
  adj_[0].resize(num_v), adj_[1].resize(num_v);
  for (int i = 0; i < es.size(); ++i) {
    V v = es[i].first, u = to(es[i].second);
    adj_[0][v].push_back(u);
    adj_[1][u].push_back(v);
  }

  for (int i = 0; i < 2; ++i) {
    for (int v = 0; v < num_v; ++v) {
      index_t tmp;
      tmp.spt_v = NULL;
      tmp.spt_d = NULL;
      tmp.spt_l = 0;
      idx_[i].push_back(tmp);
    }
  }

  //
  // Order vertices by decreasing order of degree
  //
  std::vector<V> &inv = ord_;
  std::vector<std::vector<V>> relabelled_adj[2];
  {
    // Order
    std::vector<std::pair<float, V>> deg(num_v);
    for (V v = 0; v < num_v; ++v) {
      // We add a random value here to diffuse nearby vertices
      deg[v] = std::make_pair(
          adj_[0][v].size() + adj_[1][v].size() + float(rand()) / RAND_MAX, v);
    }
    std::sort(deg.rbegin(), deg.rend());
    for (int i = 0; i < num_v; ++i) {
      V v = deg[i].second;
      if (!adj_[v].empty() || !adj_[1][v].empty()) {
        // We ignore isolated vertices here (to decide the ordering later)
        inv.push_back(v);
      }
    }

    // Relabel the vertex IDs
    relabelled_adj[0].resize(num_v), relabelled_adj[1].resize(num_v);
    std::vector<V> rank(num_v);
    for (int i = 0; i < num_v; ++i) rank[deg[i].second] = i;
    for (V from = 0; from < num_v; ++from) {
      for (size_t i = 0; i < adj_[0][from].size(); ++i) {
        V to_v = adj_[0][from][i];
        relabelled_adj[0][rank[from]].push_back(rank[to_v]);
        relabelled_adj[1][rank[to_v]].push_back(rank[from]);
      }
    }
  }

  //
  // TODO: Bit-parallel labeling
  //

  //
  // Pruned labeling
  //
  int max_threads = 1;
  for (int x = 0; x < 2; ++x) {
    const auto &adj = relabelled_adj[x];
    auto &idx = idx_[x];

    std::vector<bool> tmp_usd(num_v, false);

    // Sentinel (num_v, INF8) is added to all the vertices
    std::vector<std::pair<std::vector<V>, std::vector<W>>> tmp_idx(
        num_v,
        std::make_pair(std::vector<V>(1, INF32), std::vector<W>(1, INF8)));

    // std::vector<bool> vis(num_v);
    std::vector<uint8_t> vis(num_v / 8 + 1);
    std::vector<V> que(num_v);
    std::vector<W> dst_r(num_v + 1, INF8);

    parallel_vector<int> pdiff_nxt_que(max_threads, num_v);

    // Pruned BFS
    for (V r = 0; r < num_v; ++r) {
      if (tmp_usd[r] || adj[r].empty()) continue;
      // index_t &idx_r = index_[inv[r]];
      const auto &tmp_idx_r = tmp_idx[r];
      for (size_t i = 0; i < tmp_idx_r.first.size() - 1; ++i) {
        dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
      }

      int que_t0 = 0, que_t1 = 0, que_h = 0;
      que[que_h++] = r;
      // vis[r] = true;
      vis[r / 8] |= 1 << (r % 8);
      que_t1 = que_h;

      for (W dist = 0; que_t0 < que_h; ++dist) {
        for (int que_i = que_t0; que_i < que_t1; ++que_i) {
          V v = que[que_i];
          auto &tmp_idx_v = tmp_idx[v];
          assert(v < (V)inv.size());

          // TODO: Prefetch

          // Prune?
          if (tmp_usd[v]) continue;

          // TODO: Bit-parallel

          for (size_t i = 0; i < tmp_idx_v.first.size() - 1; ++i) {
            int w = tmp_idx_v.first[i];
            W tmp_dist = tmp_idx_v.second[i] + dst_r[w];
            if (tmp_dist <= dist) goto pruned;
          }

          // Traverse
          tmp_idx_v.first.back() = r;
          tmp_idx_v.second.back() = dist;
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
              if (__sync_bool_compare_and_swap(&vis[w / 8], m,
                                               m | (1 << (w % 8)))) {
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
      tmp_usd[r] = true;
    }

    for (int v = 0; v < (int)inv.size(); ++v) {
      int k = tmp_idx[v].first.size();
      idx[inv[v]].spt_v = (uint32_t *)memalign(64, k * sizeof(uint32_t));
      idx[inv[v]].spt_d = (uint8_t *)memalign(64, k * sizeof(uint8_t));
      idx[inv[v]].spt_l = k;
      if (!idx[inv[v]].spt_v || !idx[inv[v]].spt_d) {
        free_all();
        return;
      }
      for (int i = 0; i < k; ++i) idx[inv[v]].spt_v[i] = tmp_idx[v].first[i];
      for (int i = 0; i < k; ++i) idx[inv[v]].spt_d[i] = tmp_idx[v].second[i];
      tmp_idx[v].first.clear();
      tmp_idx[v].second.clear();
    }
    for (V i = 0; i < num_v; ++i)
      if (adj_[0][i].empty() && adj_[1][i].empty()) {
        assert(idx[i].spt_v == NULL);
        idx[i].spt_v = (uint32_t *)memalign(64, 1 * sizeof(uint32_t));
        idx[i].spt_d = (uint8_t *)memalign(64, 1 * sizeof(uint8_t));
        idx[i].spt_v[0] = INF32;
      }
  }
}

template<size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>
::query_distance(const G &g, V v_from, V v_to) {
  if (v_from == v_to) return 0;
  if (v_from >= num_v_ || v_to >= num_v_) return kInfW;

  const index_t &idx_from = idx_[1][v_from];
  const index_t &idx_to = idx_[0][v_to];
  W dist = INF8;

  // TODO: Bit-parallel

  for (int i1 = 0, i2 = 0;;) {
    V v1 = idx_from.spt_v[i1], v2 = idx_to.spt_v[i2];
    if (v1 == v2) {
      if (v1 == INF32) break;  // Sentinel
      W tmp_dist = idx_from.spt_d[i1] + idx_to.spt_d[i2];
      if (tmp_dist < dist) dist = tmp_dist;
      ++i1;
      ++i2;
    } else {
      i1 += v1 < v2 ? 1 : 0;
      i2 += v1 > v2 ? 1 : 0;
    }
  }

  if (dist >= INF8 - 2) dist = kInfW;
  return dist;
}

template<size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>
::add_edge(const G &g, V v_from, const E &e) {
  V v_to = to(e);
  V new_num_v = g.num_vertices();
  for (V v = num_v_; v < new_num_v; ++v) {
    // Now we cannot add new vertex.
    // TODO: Add index
    puts("sorry (add_edge)");
    assert(false);
  }
  num_v_ = new_num_v;

  adj_[0].resize(num_v_);
  adj_[1].resize(num_v_);
  adj_[0][v_from].push_back(v_to);
  adj_[1][v_to].push_back(v_from);

  for (int i = 0; i < 2; ++i) {
    V x = i == 0 ? v_from : v_to;
    if (adj_[0][x].size() + adj_[1][x].size() == 1) {
      // New comer
      // Now we cannot add new vertex.
      puts("sorry (add_edge)");
      assert(false);
    }
  }

  // TODO: Update bit-parallel labels
  for (int k = 0; k < kNumBitParallelRoots; ++k) {
    // PartialBPBFS(k, v_from, v_to);
  }

  //
  // Pruned BFS
  //
  {
    // v_fromに来る頂点を列挙する
    const index_t &index_from = idx_[0][v_from];
    for (int i = 0; index_from.spt_v[i] != INF32; ++i) {
      V vi = index_from.spt_v[i];
      W di = index_from.spt_d[i];
      partial_bfs(vi, v_to, di + 1, 0);
    }
  }
  {
    const index_t &index_to = idx_[1][v_to];
    for (int i = 0; index_to.spt_v[i] != INF32; ++i) {
      V vi = index_to.spt_v[i];
      W di = index_to.spt_d[i];
      partial_bfs(vi, v_from, di + 1, 1);
    }
  }
}

template<size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>
::partial_bfs(V bfs_i, V sv, W sd, int x) {
  static std::vector<std::pair<V, W>> que;
  static std::vector<bool> vis;
  // static std::vector<W> root_label;
  if ((int)que.size() < num_v_) {
    que.resize(num_v_ * 2);
    vis.resize(num_v_ * 2);
    // root_label.resize(num_v_ * 2);
  }

  // auto &index = idx_[x ^ 1];
  // V root = ord_[bfs_i];
  // index_t &idx_r = index[root];
  // for (int i = 0; idx_r.spt_v[i] != INF32; ++i) {
  //   root_label[idx_r.spt_v[i]] = idx_r.spt_d[i];
  // }

  int que_h = 0, que_t = 0;
  // queue<pair<int, int> > que;
  que[que_t++] = std::make_pair(sv, sd);
  vis[sv] = true;

  while (que_h < que_t) {
    V v = que[que_h].first;
    W d = que[que_h].second;
    ++que_h;

    // Puruning test & new label
    {
      index_t &idx_v = idx_[x][v];

      // TODO: bit-parallel

      // Case 1: |bfs_i| is already in |label[v]|
      int i = 0;
      for (; idx_v.spt_v[i] <= bfs_i; ++i) {
        V vi = idx_v.spt_v[i];
        W di = idx_v.spt_d[i];

        if (vi == bfs_i) {
          if (di <= d) {
            goto prune;
          } else {
            idx_v.spt_d[i] = d;
            goto traverse;
          }
        }
        // else if (root_label[vi] != -1 && root_label[vi] + di <= d) {
        //   goto prune;
        // }
      }

      // Case 2: |bfs_i| is not present in |label[v]|
      {
        // Expand label
        if (idx_v.spt_v[idx_v.spt_l - 1] == INF32) idx_v.Expand();
        int j;
        for (j = i; idx_v.spt_v[j] != INF32; ++j)
          ;
        for (; j >= i; --j) {
          idx_v.spt_v[j + 1] = idx_v.spt_v[j];
          idx_v.spt_d[j + 1] = idx_v.spt_d[j];
        }
      }
      idx_v.spt_v[i] = bfs_i;
      idx_v.spt_d[i] = d;
    }

  traverse:
    ;
    for (int i = 0; i < adj_[x][v].size(); ++i) {
      int w = adj_[x][v][i];
      if (!vis[w]) {
        que[que_t++] = std::make_pair(w, d + 1);
        vis[w] = true;
      }
    }

  prune:
    ;
  }

  // Reset
  for (int i = 0; i < que_t; ++i) vis[que[i].first] = false;
  // for (int i = 0; idx_r.spt_v[i] != INF32; ++i) {
  //   root_label[idx_r.spt_v[i]] = -1;
  // }
}
}  // namespace agl
