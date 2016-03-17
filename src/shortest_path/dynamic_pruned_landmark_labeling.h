#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"
#include <set>
#include <queue>

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

  static const W W_INF = 100;
  struct index_t {
    std::vector<std::pair<V, W>> spt;
    size_t size() const { return spt.size(); }
  };

  std::vector<index_t> idx[2];

  std::vector<std::vector<V>> adj[2];
  std::vector<V> inv;

 private:
  void pruned_BFS(const G &g, V v_from, int direction);
};

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::construct(
    const G &g) {
  V num_v = g.num_vertices();
  assert(num_v >= 3);
  adj[0].resize(num_v), adj[1].resize(num_v);
  idx[0].resize(num_v), idx[1].resize(num_v);

  for (const auto &p : g.edge_list()) {
    adj[0][p.first].push_back(p.second);
    adj[1][p.second].push_back(p.first);
  }

  inv.resize(num_v);
  {
    std::vector<std::pair<double, V>> decreasing_order(num_v);
    for (V i = 0; i < num_v; ++i) {
      double deg = adj[0][i].size() + adj[1][i].size();
      double t = (double)agl::random(num_v) / num_v;
      decreasing_order[i] = {deg + t, i};
    }
    std::sort(decreasing_order.rbegin(), decreasing_order.rend());
    for (int i = 0; i < num_v; ++i) inv[i] = decreasing_order[i].second;
  }

  // Pruned labelling
  for (int direction = 0; direction < 2; ++direction) {
    // decreasing order
    for (const auto &r : inv) pruned_BFS(g, r, direction);
  }
  for (int direction = 0; direction < 2; ++direction) {
    for (V v = 0; v < num_v; ++v)
      sort(idx[direction][v].spt.begin(), idx[direction][v].spt.end());
  }

  // DEBUG
  for (int direction = 0; direction < 2; ++direction) {
    std::cerr << direction << std::endl;
    for (int v = 0; v < num_v; ++v) {
      for (const auto &p : idx[direction][v].spt) {
        std::cerr << "[" << p.first << "," << p.second << "] ";
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::pruned_BFS(
    const G &g, V v_from, int direction) {
  V num_v = g.num_vertices();

  std::queue<V> que;
  que.push(v_from);

  std::vector<W> P(num_v, W_INF);
  P[v_from] = 0;

  std::vector<std::pair<V, W>> tmp_idx;

  while (!que.empty()) {
    V u = que.front();
    que.pop();
    if (query_distance(g, v_from, u) <= P[u]) continue;
    tmp_idx.push_back({u, P[u]});

    for (const auto &w : adj[direction][u]) {
      if (P[w] < W_INF) continue;
      P[w] = P[u] + 1;
      que.push(w);
    }
  }

  int another = direction ^ 1;
  for (const auto &p : tmp_idx)
    idx[another][p.first].spt.push_back({v_from, p.second});
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to) {
  if (v_from >= g.num_vertices() || v_to >= g.num_vertices())
    return v_from == v_to ? 0 : W_INF;

  W d = W_INF;
  const index_t &idx_from = idx[0][v_from];
  const index_t &idx_to = idx[1][v_to];

  // TODO bit-parallel

  for (int i1 = 0, i2 = 0; i1 < idx_from.size() && i2 < idx_to.size();) {
    V v1 = idx_from.spt[i1].first;
    V v2 = idx_to.spt[i2].first;
    if (v1 == v2) {
      W td = idx_from.spt[i1].second + idx_to.spt[i2].second;
      if (td < d) d = td;
      i1++;
      i2++;
    } else {
      i1 += v1 < v2 ? 1 : 0;
      i2 += v1 > v2 ? 1 : 0;
    }
  }
  if (d >= W_INF - 2) d = W_INF;
  return d;
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(
    const G &g, V v_from, const E &e) {
  puts("new edge");
}
}  // namespace agl
