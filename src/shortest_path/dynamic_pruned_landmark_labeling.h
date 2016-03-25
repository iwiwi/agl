#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"
#include <set>
#include <queue>
#include <bitset>

namespace agl {
template <size_t kNumBitParallelRoots = 16>
class dynamic_pruned_landmark_labeling
    : public dynamic_graph_index_interface<G>,
      public distance_query_interface<G> {
 public:
  virtual ~dynamic_pruned_landmark_labeling() {}
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

  const uint8_t D_INF = 100;
  struct index_t {
    uint8_t bpspt_d[kNumBitParallelRoots];
    uint64_t bpspt_s[kNumBitParallelRoots][2];
    std::vector<std::pair<V, uint8_t>> spt_p;
    void update(V v, uint8_t d) {
      if (spt_p.empty() || spt_p.rbegin()->first < v) {
        spt_p.emplace_back(v, d);
        return;
      }
      auto it = std::lower_bound(spt_p.begin(), spt_p.end(),
                                 std::make_pair(v, (uint8_t)0));
      if (it->first == v) {
        if (it->second > d) it->second = d;
        return;
      }
      if (it->first > v) {
        spt_p.insert(it, std::make_pair(v, d));
        return;
      }
      assert(false);
    }
    size_t size() const { return spt_p.size(); }
  };

  std::vector<index_t> idx[2];
  std::vector<std::vector<V>> adj[2];
  std::vector<V> rank;
  std::vector<V> inv;

  size_t total_label_num() {
    size_t sum = 0;
    for (int i = 0; i < 2; ++i)
      for (const index_t &j : idx[i]) sum += j.size();
    return sum;
  }

  // private:
  void init_graph(const G &g);
  void bfs_together(V root, const std::vector<bool> &used);
  void pruned_bfs(V root, int direction, const std::vector<bool> &used);
  void resume_pbfs(V v_from, V v_to, uint8_t d_ft, int direction);
  uint8_t query_distance_(V v_from, V v_to, int direction);
  void bit_parallel_bfs(std::vector<bool> &used);
  void partial_bp_bfs(int bp_i, V v_from, int direction);
};

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::construct(
    const G &g) {
  V num_v = g.num_vertices();
  assert(num_v >= 3);

  // Initialize
  init_graph(g);

  // Bit-Parallel Labeling
  std::vector<bool> used(num_v, false);
  bit_parallel_bfs(used);

  // Pruned labelling
  for (V root = 0; root < num_v; ++root) {
    if (used[root]) continue;
    pruned_bfs(root, 0, used);
    pruned_bfs(root, 1, used);
    used[root] = true;
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::bit_parallel_bfs(
    std::vector<bool> &used) {
  const V num_v = rank.size();
  V root = 0;
  for (int bp_i = 0; bp_i < kNumBitParallelRoots; ++bp_i) {
    // Select Root
    while (root < num_v && used[root]) root++;
    if (root == num_v) {
      for (V v = 0; v < num_v; ++v) {
        idx[0][v].bpspt_d[bp_i] = D_INF;
        idx[1][v].bpspt_d[bp_i] = D_INF;
      }
      continue;
    }
    used[root] = true;

    // Select Roots
    std::vector<V> selected_roots;
    const std::vector<V> &adjo = adj[0][root];
    const std::vector<V> &adji = adj[1][root];
    for (int o = 0, i = 0; o < adjo.size() && i < adji.size();) {
      V vo = adjo[o], vi = adji[i];
      if (vo == vi) {
        if (!used[vo]) {
          used[vo] = true;
          selected_roots.push_back(vo);
          if (selected_roots.size() == 64) break;
        }
        o++;
        i++;
      } else {
        if (vo > vi) i++;
        if (vo < vi) o++;
      }
    }

    for (int direction = 0; direction < 2; ++direction) {
      int another = direction ^ 1;

      std::vector<uint8_t> tmp_d(num_v, D_INF);
      std::vector<std::pair<uint64_t, uint64_t>> tmp_s(num_v, {0, 0});
      std::queue<std::pair<V, uint8_t>> que;

      que.emplace(root, 0);
      tmp_d[root] = 0;

      for (size_t i = 0; i < selected_roots.size(); ++i) {
        V v = selected_roots[i];
        que.emplace(v, 1);
        tmp_d[v] = 1;
        tmp_s[v].first = 1ULL << i;
      }

      for (uint8_t d = 0; !que.empty(); ++d) {
        std::vector<std::pair<V, V>> sibling_es;
        std::vector<std::pair<V, V>> child_es;
        while (!que.empty() && que.front().second == d) {
          V v = que.front().first;
          que.pop();
          for (V tv : adj[direction][v]) {
            uint8_t td = d + 1;
            if (d > tmp_d[tv]) continue;
            if (d == tmp_d[tv] && v < tv) sibling_es.emplace_back(v, tv);
            if (d < tmp_d[tv]) child_es.emplace_back(v, tv);
            if (tmp_d[tv] == D_INF) {
              tmp_d[tv] = td;
              que.emplace(tv, td);
            }
          }
        }

        for (auto &p : sibling_es) {
          V v, w;
          std::tie(v, w) = p;
          tmp_s[v].second |= tmp_s[w].first;
          tmp_s[w].second |= tmp_s[v].first;
        }
        for (auto &p : child_es) {
          V v, c;
          std::tie(v, c) = p;
          tmp_s[c].first |= tmp_s[v].first;
          tmp_s[c].second |= tmp_s[v].second;
        }
      }

      for (V v = 0; v < num_v; ++v) {
        idx[another][v].bpspt_d[bp_i] = tmp_d[v];
        idx[another][v].bpspt_s[bp_i][0] = tmp_s[v].first;
        idx[another][v].bpspt_s[bp_i][1] = tmp_s[v].second & ~tmp_s[v].first;
      }
    }
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::init_graph(
    const G &g) {
  V num_v = g.num_vertices();
  rank.resize(num_v);
  inv.resize(num_v);
  {
    std::vector<V> deg(num_v, 0);
    for (const auto &p : g.edge_list())
      if (p.first != p.second) deg[p.first]++, deg[p.second]++;

    std::vector<std::pair<double, V>> sorting_v(num_v);
    for (V v = 0; v < num_v; ++v) {
      double t = (double)agl::random(num_v) / num_v;
      sorting_v[v] = std::make_pair(deg[v] + t, v);
    }
    std::sort(sorting_v.rbegin(), sorting_v.rend());
    for (int i = 0; i < num_v; ++i) inv[i] = sorting_v[i].second;
    for (int i = 0; i < num_v; ++i) rank[inv[i]] = i;
  }

  idx[0].resize(num_v), idx[1].resize(num_v);
  adj[0].resize(num_v), adj[1].resize(num_v);
  for (const auto &p : g.edge_list()) {
    V v = p.first, u = p.second;
    if (v == u) continue;
    adj[0][rank[v]].push_back(rank[u]);
    adj[1][rank[u]].push_back(rank[v]);
  }
  for (V v = 0; v < num_v; ++v) {
    std::sort(adj[0][v].begin(), adj[0][v].end());
    std::sort(adj[1][v].begin(), adj[1][v].end());
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::pruned_bfs(
    V root, int direction, const std::vector<bool> &used) {
  int another = direction ^ 1;
  V num_v = rank.size();
  std::queue<V> que;
  que.push(root);
  std::vector<uint8_t> P(num_v, D_INF);
  P[root] = 0;

  while (!que.empty()) {
    V u = que.front();
    que.pop();
    if (used[u]) continue;
    if (u != root && query_distance_(root, u, direction) <= P[u]) continue;
    idx[another][u].update(root, P[u]);

    for (const auto &w : adj[direction][u]) {
      if (P[w] < D_INF) continue;
      P[w] = P[u] + 1;
      que.push(w);
    }
  }
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to) {
  v_from = rank[v_from], v_to = rank[v_to];
  return query_distance_(v_from, v_to, 0);
}

template <size_t kNumBitParallelRoots>
uint8_t dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance_(
    V v_from, V v_to, int direction) {
  int another = direction ^ 1;
  assert(v_from < adj[0].size() && v_to < adj[0].size());
  if (v_from == v_to) return 0;

  uint8_t d = D_INF;
  const index_t &idx_from = idx[direction][v_from];
  const index_t &idx_to = idx[another][v_to];

  for (int bp_i = 0; bp_i < kNumBitParallelRoots; ++bp_i) {
    uint8_t td = idx_from.bpspt_d[bp_i] + idx_to.bpspt_d[bp_i];
    if (td - 2 > d) continue;
    auto &from_masks = idx_from.bpspt_s[bp_i];
    auto &to_masks = idx_to.bpspt_s[bp_i];
    if (from_masks[0] & to_masks[0]) {
      td -= 2;
    } else if ((from_masks[1] & to_masks[0]) | (from_masks[0] & to_masks[1])) {
      td -= 1;
    }
    if (td < d) d = td;
  }

  for (auto &p : idx_from.spt_p)
    if (p.first == v_to && d > p.second) d = p.second;
  for (auto &p : idx_to.spt_p)
    if (p.first == v_from && d > p.second) d = p.second;

  for (auto i1 = idx_from.spt_p.begin(), i2 = idx_to.spt_p.begin();
       i1 != idx_from.spt_p.end() && i2 != idx_to.spt_p.end();) {
    V v1 = i1->first;
    V v2 = i2->first;
    if (v1 == v2) {
      uint8_t td = i1->second + i2->second;
      if (td < d) d = td;
      i1++;
      i2++;
    } else {
      if (v1 < v2) i1++;
      if (v1 > v2) i2++;
    }
  }
  if (d >= D_INF - 2) d = D_INF;
  return d;
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::resume_pbfs(
    V neighbor_a, V v_b, uint8_t d_nab, int direction) {
  int another = direction ^ 1;

  std::queue<std::pair<V, uint8_t>> que;
  que.emplace(v_b, d_nab);
  while (!que.empty()) {
    V v;
    uint8_t d;
    std::tie(v, d) = que.front();
    que.pop();
    if (query_distance_(neighbor_a, v, direction) <= d) continue;
    idx[another][v].update(neighbor_a, d);
    for (const auto &w : adj[direction][v]) que.emplace(w, d + 1);
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::partial_bp_bfs(
    int bp_i, V v_from, int direction) {
  V another = direction ^ 1;
  const index_t &idx_from = idx[another][v_from];
  const uint8_t base_d = idx_from.bpspt_d[bp_i];
  if (base_d == D_INF) return;

  std::queue<std::pair<V, uint8_t>> que;
  que.emplace(v_from, base_d);
  for (uint8_t d = base_d; !que.empty(); ++d) {
    std::vector<V> tmp;
    while (!que.empty() && que.front().second == d) {
      V v = que.front().first;
      const index_t &idx_v = idx[another][v];
      que.pop();
      tmp.push_back(v);

      for (V tv : adj[direction][v]) {
        index_t &idx_tv = idx[another][tv];
        if (d == idx_tv.bpspt_d[bp_i])
          idx_tv.bpspt_s[bp_i][1] |= idx_v.bpspt_s[bp_i][0];
        if (d + 1 < idx_tv.bpspt_d[bp_i]) {
          idx_tv.bpspt_s[bp_i][1] =
              (idx_tv.bpspt_d[bp_i] == d + 2 ? idx_tv.bpspt_s[bp_i][0] : 0);

          idx_tv.bpspt_s[bp_i][0] = 0;
          idx_tv.bpspt_d[bp_i] = d + 1;
          que.emplace(tv, d + 1);
        }
      }
    }

    for (V v : tmp) {
      const index_t &idx_v = idx[another][v];
      for (V tv : adj[direction][v]) {
        index_t &idx_tv = idx[another][tv];
        if (idx_tv.bpspt_d[bp_i] == d + 1) {
          // Set propagation (2)
          idx_tv.bpspt_s[bp_i][0] |= idx_v.bpspt_s[bp_i][0];
          idx_tv.bpspt_s[bp_i][1] |= idx_v.bpspt_s[bp_i][1];
        }
      }
    }
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(
    const G &g, V v_a, const E &e) {
  V v_b = to(e);
  v_a = rank[v_a], v_b = rank[v_b];
  assert(v_a < g.num_vertices() && v_b < g.num_vertices());
  if (query_distance_(v_a, v_b, 0) <= 1) return;

  adj[0][v_a].push_back(v_b);
  adj[1][v_b].push_back(v_a);
  std::sort(adj[0][v_a].begin(), adj[0][v_a].end());
  std::sort(adj[1][v_b].begin(), adj[1][v_b].end());

  for (int bp_i = 0; bp_i < kNumBitParallelRoots; ++bp_i) {
    partial_bp_bfs(bp_i, v_a, 0);
    partial_bp_bfs(bp_i, v_b, 1);
  }

  std::vector<std::tuple<V, uint8_t, int>> tmp;
  for (const auto &p : idx[1][v_a].spt_p)
    tmp.emplace_back(p.first, p.second, 0);
  for (const auto &p : idx[0][v_b].spt_p)
    tmp.emplace_back(p.first, p.second, 1);
  std::sort(tmp.begin(), tmp.end());

  for (const auto &q : tmp) {
    V r, t;
    uint8_t d;
    std::tie(r, d, t) = q;
    if (t == 0) {
      resume_pbfs(r, v_b, d + 1, 0);
    }
    if (t == 1) {
      resume_pbfs(r, v_a, d + 1, 1);
    }
  }
}
}  // namespace agl
