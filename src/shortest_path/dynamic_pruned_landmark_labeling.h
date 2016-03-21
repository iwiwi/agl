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
  // private:
  static const W W_INF = 100;
  struct index_t {
    std::map<V, W> spm;
    void update(V v, W d) {
      if (spm.count(v) && spm[v] <= d) return;
      spm[v] = d;
    }
  };

  std::vector<index_t> idx[2];
  std::vector<std::vector<V>> adj[2];
  std::vector<V> rank;
  std::vector<V> inv;

 private:
  void pruned_bfs(V root, int direction);
  void resume_pbfs(V v_from, V v_to, W d_ft, int direction);
  W query_distance(V v_from, V v_to, int direction);
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

  for (V v = 0; v < num_v; ++v) {
    idx[0][v].update(v, 0);
    idx[1][v].update(v, 0);
  }

  rank.resize(num_v);
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
    for (int i = 0; i < num_v; ++i) rank[inv[i]] = i;
  }

  // Pruned labelling
  for (const auto &r : inv) {
    pruned_bfs(r, 0);
    pruned_bfs(r, 1);
  }

  // // DEBUG
  // for (int dir = 0; dir < 2; ++dir) {
  //   for (int v = 0; v < num_v; ++v) {
  //     std::cerr << v << " ";
  //     for (const auto &p : idx[dir][v].spm)
  //       std::cerr << "[" << p.first << "," << p.second << "] ";
  //     std::cerr << std::endl;
  //   }
  // }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::pruned_bfs(
    V root, int direction) {
  V num_v = adj[0].size();
  std::queue<V> que;
  que.push(root);
  std::vector<W> P(num_v, W_INF);
  P[root] = 0;

  std::vector<std::pair<V, W>> tmp_idx;

  while (!que.empty()) {
    V u = que.front();
    que.pop();
    if (u != root && query_distance(root, u, direction) <= P[u]) continue;
    tmp_idx.emplace_back(u, P[u]);

    for (const auto &w : adj[direction][u]) {
      if (P[w] < W_INF) continue;
      P[w] = P[u] + 1;
      que.push(w);
    }
  }

  int another = direction ^ 1;
  for (const auto &p : tmp_idx) idx[another][p.first].update(root, p.second);
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to) {
  return query_distance(v_from, v_to, 0);
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    V v_from, V v_to, int direction) {
  int another = direction ^ 1;
  assert(v_from < adj[0].size() && v_to < adj[0].size());
  if (v_from == v_to) return 0;

  W d = W_INF;
  const index_t &idx_from = idx[direction][v_from];
  const index_t &idx_to = idx[another][v_to];

  // TODO bit-parallel

  for (auto i1 = idx_from.spm.begin(), i2 = idx_to.spm.begin();
       i1 != idx_from.spm.end() && i2 != idx_to.spm.end();) {
    V v1 = i1->first;
    V v2 = i2->first;
    if (v1 == v2) {
      W td = i1->second + i2->second;
      if (td < d) d = td;
      i1++;
      i2++;
    } else {
      if (v1 < v2) i1++;
      if (v1 > v2) i2++;
    }
  }
  if (d >= W_INF - 2) d = W_INF;
  return d;
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::resume_pbfs(
    V neighbor_a, V v_b, W d_nab, int direction) {
  int another = direction ^ 1;

  std::queue<std::pair<V, W>> que;
  que.emplace(v_b, d_nab);
  while (!que.empty()) {
    V v;
    W d;
    std::tie(v, d) = que.front();
    que.pop();
    if (query_distance(neighbor_a, v, direction) <= d) continue;
    idx[another][v].update(neighbor_a, d);
    for (const auto &w : adj[direction][v]) que.emplace(w, d + 1);
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(
    const G &g, V v_a, const E &e) {
  V v_b = to(e);
  assert(v_a < g.num_vertices() && v_b < g.num_vertices());
  if (query_distance(g, v_a, v_b) <= 1) return;
  adj[0][v_a].push_back(v_b);
  adj[1][v_b].push_back(v_a);

  std::vector<std::tuple<V, W, int>> tmp;
  for (const auto &p : idx[1][v_a].spm) {
    V r = rank[p.first];
    V d = p.second;
    tmp.emplace_back(r, d, 0);
  }
  for (const auto &p : idx[0][v_b].spm) {
    V r = rank[p.first];
    V d = p.second;
    tmp.emplace_back(r, d, 1);
  }

  std::sort(tmp.begin(), tmp.end());

  for (const auto &q : tmp) {
    V r, t;
    W d;
    std::tie(r, d, t) = q;
    if (t == 0) {
      V neighbor_a = inv[r];
      resume_pbfs(neighbor_a, v_b, d + 1, 0);
    }
    if (t == 1) {
      V neighbor_b = inv[r];
      resume_pbfs(neighbor_b, v_a, d + 1, 1);
    }
  }
}
}  // namespace agl
