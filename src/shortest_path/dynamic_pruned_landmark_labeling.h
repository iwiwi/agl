#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"
#include <malloc.h>
#include <bitset>

namespace agl {
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
  const W INF8 = 100;  // For unreachable pairs
  const uint32_t INF32 = std::numeric_limits<int32_t>::max();  // For sentinel

  struct index_t {
    W bpspt_d[kNumBitParallelRoots];
    uint64_t bpspt_s[kNumBitParallelRoots][2];
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
  } __attribute__((aligned(64)));  // Aligned for cache lines

  V num_v_;
  std::vector<std::vector<V>> adj_[2];
  std::vector<index_t> idx_[2];
  std::vector<V> ord_;
  // bp_inv_[x][i].first=ordered_root => second=i_bpspt
  std::vector<std::pair<V, size_t>> bp_inv_[2];

  // ord_negighbors_[v].first:= vに隣接する頂点,
  // ord_negighbors_[v].second:=v->firstならtrue
  std::vector<std::vector<std::pair<V, bool>>> ord_negighbors_;

  void free_all() {
    for (int i = 0; i < 2; ++i) {
      for (V v = 0; v < num_v_; ++v) {
        free(idx_[i][v].spt_v);
        idx_[i][v].spt_v = NULL;
        free(idx_[i][v].spt_d);
        idx_[i][v].spt_d = NULL;
      }
      idx_[i].clear();
    }
    num_v_ = 0;
  }

  void print_index() {
    for (int x = 0; x < 2; ++x) {
      std::cerr << "Index No. " << x << std::endl;
      for (V v = 0; v < num_v_; ++v) {
        const index_t &idx = idx_[x][v];
        for (int i = 0; idx.spt_v[i] != INF32; ++i) {
          W d = idx.spt_d[i];
          V ordered = idx.spt_v[i];
          std::cerr << v << "->" << ord_[ordered] << " " << d << " " << ordered
                    << std::endl;
        }
      }
    }
  }
  void partial_bfs(V bfs_i, V sv, W sd, int x);
};

template<size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>
::construct(const G &g) {
  free_all();
  V &num_v = num_v_;
  num_v = g.num_vertices();
  assert(num_v >= 3);
  const std::vector<std::pair<V, E>> &es = g.edge_list();
  adj_[0].resize(num_v), adj_[1].resize(num_v);
  ord_negighbors_.resize(num_v);
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
      if (!adj_[0][v].empty() || !adj_[1][v].empty()) {
        // We ignore isolated vertices here (to decide the ordering later)
        inv.push_back(v);
      }
    }

    // Relabel the vertex IDs
    relabelled_adj[0].resize(num_v), relabelled_adj[1].resize(num_v);
    std::vector<V> rank(num_v);
    for (int i = 0; i < num_v; ++i) rank[deg[i].second] = i;
    for (V from = 0; from < num_v; ++from) {
      for (V v_to : adj_[0][from]) {
        relabelled_adj[0][rank[from]].push_back(rank[v_to]);
        relabelled_adj[1][rank[v_to]].push_back(rank[from]);
        ord_negighbors_[rank[from]].push_back(std::make_pair(rank[v_to], true));
        ord_negighbors_[rank[v_to]].push_back(
            std::make_pair(rank[from], false));
      }
    }
  }

  //
  // Bit-parallel labeling
  //

  for (int x = 0; x < 2; ++x) {
    std::vector<std::vector<V>> &adj = relabelled_adj[x];
    std::vector<index_t> &index = idx_[x];
    std::vector<std::pair<V, size_t>> &bp_inv = bp_inv_[x];

    size_t num_e = g.num_edges();
    std::vector<bool> used(num_v, false);
    std::vector<W> tmp_d(num_v);  // Distance from root
    std::vector<std::pair<uint64_t, uint64_t>> tmp_s(num_v);
    std::vector<V> que(num_v);
    std::vector<std::pair<V, V>> sibling_es(num_e);
    std::vector<std::pair<V, V>> child_es(num_e);

    V root = 0;
    for (size_t i_bpspt = 0; i_bpspt < kNumBitParallelRoots; ++i_bpspt) {
      while (root < num_v && used[root]) {
        ++root;
      }
      if (root == num_v) {
        for (V v = 0; v < num_v; ++v) index[v].bpspt_d[i_bpspt] = INF8;
        continue;
      }
      bp_inv.push_back(std::make_pair(root, i_bpspt));

      used[root] = true;

      fill(tmp_d.begin(), tmp_d.end(), INF8);
      fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));

      int que_t0 = 0, que_t1 = 0, que_h = 0;
      que[que_h++] = root;
      tmp_d[root] = 0;
      que_t1 = que_h;

      // Start from root
      int bit_pos = 0;
      sort(ord_negighbors_[root].begin(), ord_negighbors_[root].end());
      for (std::pair<V, bool> neighbor_pair : ord_negighbors_[root]) {
        V v = neighbor_pair.first;
        if ((x == 0 && neighbor_pair.second) ||
            (x == 1 && !neighbor_pair.second)) {
          // root->v
          if (!used[v]) {
            tmp_d[v] = 1;
            tmp_s[v].first = 1ULL << bit_pos;
            used[v] = true;
            que[que_h++] = v;
          }
        }
        if (++bit_pos == 64) break;
      }

      for (W dist = 0; que_t0 < que_h; ++dist) {
        int num_sibling_es = 0, num_child_es = 0;

        for (int que_i = que_t0; que_i < que_t1; ++que_i) {
          V v = que[que_i];

          for (V tv : adj[v]) {
            if (dist == tmp_d[tv]) {
              if (v < tv) {  // to avoid duplicate pair(v,tv)
                // dist(root->v) = dist(root->tv)
                sibling_es[num_sibling_es].first = v;
                sibling_es[num_sibling_es].second = tv;
                ++num_sibling_es;
              }
            } else if (dist < tmp_d[tv]) {
              if (tmp_d[tv] == INF8) {
                que[que_h++] = tv;
                tmp_d[tv] = dist + 1;
              }
              // dist(root->v) + 1 = dist(root->tv)
              // and root->v->tv
              child_es[num_child_es].first = v;
              child_es[num_child_es].second = tv;
              ++num_child_es;
            }
          }
        }

        for (int i = 0; i < num_sibling_es; ++i) {
          V v = sibling_es[i].first, w = sibling_es[i].second;
          tmp_s[v].second |= tmp_s[w].first;
          tmp_s[w].second |= tmp_s[v].first;
        }
        for (int i = 0; i < num_child_es; ++i) {
          V parent = child_es[i].first;
          V child = child_es[i].second;
          tmp_s[child].first |= tmp_s[parent].first;
          tmp_s[child].second |= tmp_s[parent].second;
        }

        que_t0 = que_t1;
        que_t1 = que_h;
      }

      for (V v = 0; v < (V)inv.size(); ++v) {
        index[inv[v]].bpspt_d[i_bpspt] = tmp_d[v];

        // inv[v]と, inv[v]よりrootに近い頂点のbit (rootは除く)
        index[inv[v]].bpspt_s[i_bpspt][0] = tmp_s[v].first;
        // inv[v]とrootから等距離にある点のbit(inv[v]以外)
        index[inv[v]].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
      }
      for (V v = 0; v < num_v; ++v) {
        if (adj_[x][v].empty()) {
          index[v].bpspt_d[i_bpspt] = INF8;
          index[v].bpspt_s[i_bpspt][0] = index[v].bpspt_s[i_bpspt][1] = 0;
        }
      }
    }  // finish on a root

    std::sort(bp_inv.begin(), bp_inv.end());
  }

  //
  // Pruned labeling
  //
  for (int x = 0; x < 2; ++x) {
    std::vector<bool> used(num_v, false);

    // tmp_idx[v][i].second:= Distance from |tmp_idx[v][i].first| to |v|.
    // Sentinel (num_v, INF8) is added to all the vertices
    std::vector<std::vector<std::pair<V, W>>> tmp_idx(
        num_v, std::vector<std::pair<V, W>>(1, std::make_pair(INF32, INF8)));

    std::vector<bool> vis(num_v, false);
    std::vector<V> que(num_v);

    // Pruned BFS
    for (V ordered_root = 0; ordered_root < num_v; ++ordered_root) {
      if (used[ordered_root]) continue;

      int que_t0 = 0, que_t1 = 0, que_h = 0;
      que[que_h++] = ordered_root;
      vis[ordered_root] = true;
      que_t1 = que_h;

      for (W dist = 0; que_t0 < que_h; ++dist) {
        for (int que_i = que_t0; que_i < que_t1; ++que_i) {
          V v = que[que_i];
          auto &tmp_idx_v = tmp_idx[v];
          assert(v < (V)inv.size());

          // Traverse
          tmp_idx_v.back().first = ordered_root;
          tmp_idx_v.back().second = dist;
          tmp_idx_v.push_back(std::make_pair(INF32, INF8));

          // Prune?
          if (used[v]) continue;

          //
          // Bit-parallel Prune
          //
          {
            V v_from = inv[ordered_root];
            V v_to = inv[v];
            const index_t &idx_from = idx_[x ^ 1][v_from];
            const index_t &idx_to = idx_[x][v_to];

            for (size_t i = 0, j = 0;
                 i < bp_inv_[x].size() && j < bp_inv_[x ^ 1].size();) {
              V vi = bp_inv_[x][i].first;
              V vj = bp_inv_[x ^ 1][j].first;
              if (vi == vj) {
                if (vi == num_v_) break;
                size_t i_bpspt = bp_inv_[x][i].second;
                size_t j_bpspt = bp_inv_[x ^ 1][j].second;

                W tmp_dist =
                    idx_from.bpspt_d[j_bpspt] + idx_to.bpspt_d[i_bpspt];
                if (tmp_dist - 2 < dist) {
                  if (idx_from.bpspt_s[j_bpspt][x] &
                      idx_to.bpspt_s[i_bpspt][x]) {
                    tmp_dist -= 2;
                  } else if ((idx_from.bpspt_s[j_bpspt][x] &
                              idx_to.bpspt_s[i_bpspt][x ^ 1]) |
                             (idx_from.bpspt_s[j_bpspt][x ^ 1] &
                              idx_to.bpspt_s[i_bpspt][x])) {
                    tmp_dist--;
                  }
                  if (tmp_dist < dist) {
                    goto pruned;
                  }
                }
                i++;
                j++;
              } else if (vi > vj) {
                j++;
              } else {
                i++;
              }
            }
          }

          // for (size_t i = 0; i < tmp_idx_v.first.size() - 1; ++i) {
          //   int w = tmp_idx_v.first[i];
          //   int td = tmp_idx_v.second[i] + dst_r[w];
          //   if (td <= dist) goto pruned;
          // }

          for (V w : relabelled_adj[x][v]) {
            if (!vis[w]) {
              que[que_h++] = w;
              vis[w] = true;
            }
          }
        pruned : {}
        }

        que_t0 = que_t1;
        que_t1 = que_h;
      }

      // Reset for next BFS
      for (int i = 0; i < vis.size(); ++i) vis[i] = false;
      used[ordered_root] = true;
    }  // Pruned BFS

    auto &idx = idx_[x];
    for (V v = 0; v < (V)inv.size(); ++v) {
      int k = tmp_idx[v].size();
      idx[inv[v]].spt_v = (uint32_t *)memalign(64, k * sizeof(uint32_t));
      idx[inv[v]].spt_d = (uint8_t *)memalign(64, k * sizeof(uint8_t));
      idx[inv[v]].spt_l = k;
      if (!idx[inv[v]].spt_v || !idx[inv[v]].spt_d) {
        free_all();
        return;
      }
      for (int i = 0; i < k; ++i) idx[inv[v]].spt_v[i] = tmp_idx[v][i].first;
      for (int i = 0; i < k; ++i) idx[inv[v]].spt_d[i] = tmp_idx[v][i].second;
      tmp_idx[v].clear();
    }

    // Allocate for isolated vertices
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
  if (v_from >= num_v_ || v_to >= num_v_) return INF8;

  const index_t &idx_from = idx_[1][v_from];
  const index_t &idx_to = idx_[0][v_to];
  W dist = INF8;

  for (size_t i = 0, j = 0; i < bp_inv_[0].size() && j < bp_inv_[1].size();) {
    V vi = bp_inv_[0][i].first;
    V vj = bp_inv_[1][j].first;
    if (vi == vj) {
      if (vi == num_v_) break;
      size_t i_bpspt = bp_inv_[0][i].second;
      size_t j_bpspt = bp_inv_[1][j].second;

      W tmp_dist = idx_from.bpspt_d[j_bpspt] + idx_to.bpspt_d[i_bpspt];
      if (tmp_dist - 2 < dist) {
        if (idx_from.bpspt_s[j_bpspt][0] & idx_to.bpspt_s[i_bpspt][0]) {
          tmp_dist -= 2;
        } else if ((idx_from.bpspt_s[j_bpspt][0] & idx_to.bpspt_s[i_bpspt][1]) |
                   (idx_from.bpspt_s[j_bpspt][1] &
                    idx_to.bpspt_s[i_bpspt][0])) {
          tmp_dist--;
        }
        if (tmp_dist < dist) {
          dist = tmp_dist;
        }
      }
      i++;
      j++;
    } else if (vi > vj) {
      j++;
    } else {
      i++;
    }
  }

  for (int i1 = 0, i2 = 0;;) {
    V v1 = idx_from.spt_v[i1], v2 = idx_to.spt_v[i2];
    W d1 = idx_from.spt_d[i1], d2 = idx_to.spt_d[i2];
    if (v1 == v2) {
      if (v1 == INF32) break;  // Sentinel
      W tmp_dist = d1 + d2;
      if (tmp_dist < dist) dist = tmp_dist;
      ++i1;
      ++i2;
    } else {
      i1 += v1 < v2 ? 1 : 0;
      i2 += v1 > v2 ? 1 : 0;
    }
  }

  if (dist >= INF8 - 2) dist = INF8;
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
  assert(v_from >= 0 && v_to >= 0 && v_from < num_v_ && v_to < num_v_);

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

  // bit-parallel partial BFS
  {
    const index_t &idx_from = idx_[1][v_from];
    const index_t &idx_to = idx_[0][v_to];

    for (size_t i = 0, j = 0; i < bp_inv_[0].size() && j < bp_inv_[1].size();) {
      V vi = bp_inv_[0][i].first;
      V vj = bp_inv_[1][j].first;
      if (vi == vj) {
        if (vi == num_v_) break;
        size_t i_bpspt = bp_inv_[0][i].second;

        W base_dist = idx_to.bpspt_d[i_bpspt];  // vi=vj -> v_to

        static std::vector<V> que0, que1;
        static std::vector<bool> queued;
        if ((V)queued.size() < num_v_) {
          que0.resize(num_v_ * 2);
          que1.resize(num_v_ * 2);
          queued.resize(num_v_ * 2);
        }

        int que_h0 = 0, que_h1 = 0;
        queued[que0[que_h0++] = v_to] = true;

        for (W dist = base_dist; que_h0 > 0; ++dist) {
          for (int i0 = 0; i0 < que_h0; ++i0) {
            V v = que0[i0];
            W td = dist + 1;
            index_t &idx_v = idx_[0][v];
            for (V tv : adj_[0][v]) {
              index_t &idx_tv = idx_[0][tv];
              if (dist == idx_tv.bpspt_d[i_bpspt]) {
                idx_tv.bpspt_s[i_bpspt][1] |= idx_v.bpspt_s[i_bpspt][0];
              }
              if (td < idx_tv.bpspt_d[i_bpspt]) {
                if (idx_tv.bpspt_d[i_bpspt] == dist + 1 + 1) {
                  idx_tv.bpspt_s[i_bpspt][1] = idx_tv.bpspt_s[i_bpspt][0];
                } else {
                  idx_tv.bpspt_s[i_bpspt][1] = 0;
                }

                idx_tv.bpspt_s[i_bpspt][0] = 0;
                idx_tv.bpspt_d[i_bpspt] = dist + 1;
                assert(!queued[tv]);
                queued[que1[que_h1++] = tv] = true;
              }
            }
          }

          for (int i0 = 0; i0 < que_h0; ++i0) {
            V v = que0[i0];
            for (V tv : adj_[0][v]) {
              if (idx_[0][tv].bpspt_d[i_bpspt] == dist + 1) {
                // Set propagation (2)
                idx_[0][tv].bpspt_s[i_bpspt][0] |=
                    idx_[0][v].bpspt_s[i_bpspt][0];
                idx_[0][tv].bpspt_s[i_bpspt][1] |=
                    idx_[0][v].bpspt_s[i_bpspt][1];
              }
            }
          }

          for (int i0 = 0; i0 < que_h0; ++i0) queued[que0[i0]] = false;
          que0.swap(que1);
          que_h0 = que_h1;
          que1.clear();
          que_h1 = 0;
        }

        i++;
        j++;
      } else if (vi > vj) {
        j++;
      } else {
        i++;
      }
    }
  }

  //
  // Pruned BFS
  //
  {
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
  if ((int)que.size() < num_v_) {
    que.resize(num_v_ * 2);
    vis.resize(num_v_ * 2);
  }

  int que_h = 0, que_t = 0;
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
      {
        const index_t &idx_from = idx_[x ^ 1][ord_[bfs_i]];

        for (size_t i = 0, j = 0;
             i < bp_inv_[x].size() && j < bp_inv_[x ^ 1].size();) {
          V vi = bp_inv_[x][i].first;
          V vj = bp_inv_[x ^ 1][j].first;
          if (vi == vj) {
            if (vi == num_v_) break;
            size_t i_bpspt = bp_inv_[x][i].second;
            size_t j_bpspt = bp_inv_[x ^ 1][j].second;

            W tmp_dist = idx_from.bpspt_d[j_bpspt] + idx_v.bpspt_d[i_bpspt];
            if (tmp_dist - 2 < d) {
              if (idx_from.bpspt_s[j_bpspt][x] & idx_v.bpspt_s[i_bpspt][x]) {
                tmp_dist -= 2;
              } else if ((idx_from.bpspt_s[j_bpspt][x] &
                          idx_v.bpspt_s[i_bpspt][x ^ 1]) |
                         (idx_from.bpspt_s[j_bpspt][x ^ 1] &
                          idx_v.bpspt_s[i_bpspt][x])) {
                tmp_dist--;
              }
              if (tmp_dist < d) {
                goto prune;
              }
            }
            i++;
            j++;
          } else if (vi > vj) {
            j++;
          } else {
            i++;
          }
        }
      }
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
        } else {
          const index_t &idx_k = idx_[x][ord_[vi]];
          int k = 0;
          for (; idx_k.spt_v[k] < bfs_i; ++k)
            ;
          if (idx_k.spt_v[k] == bfs_i) {
            W tmp_d = idx_k.spt_d[k] + di;
            if (tmp_d <= d) goto prune;
          }
        }
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
}

}  // namespace agl
