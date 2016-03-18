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
    void update(V v, W d) {
      auto it = std::lower_bound(
          spt.begin(), spt.end(), std::make_pair(v, -1),
          [](std::pair<int, int> lp, std::pair<int, int> rp) -> bool {
            if (lp.first == rp.first) return lp.second < rp.second;
            return lp.first < rp.first;
          });
      if (it == spt.end()) {
        spt.emplace_back(v, d);
        return;
      }
      if (it->first == v) {
        if (it->second > d) it->second = d;
        return;
      }
      assert(it->first > v);
      spt.insert(it, {v, d});
    }
  };

  std::vector<index_t> idx[2];

  std::vector<std::vector<V>> adj[2];

 private:
  void pruned_bfs(const G &g, V v_from, int direction);
  void resume_pbfs(const G &g, V v_from, V v_to, W d_ft, int direction);
  W query_distance(const G &g, V v_from, V v_to, int direction);
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

  std::vector<V> inv(num_v);
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
    for (const auto &r : inv) pruned_bfs(g, r, direction);
  }
  for (int direction = 0; direction < 2; ++direction) {
    for (V v = 0; v < num_v; ++v)
      sort(idx[direction][v].spt.begin(), idx[direction][v].spt.end());
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::pruned_bfs(
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
  return query_distance(g, v_from, v_to, 0);
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to, int direction) {
  int another = direction ^ 1;
  if (v_from >= g.num_vertices() || v_to >= g.num_vertices())
    return v_from == v_to ? 0 : W_INF;

  W d = W_INF;
  const index_t &idx_from = idx[direction][v_from];
  const index_t &idx_to = idx[another][v_to];

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
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::resume_pbfs(
    const G &g, V v_from, V v_to, W d_ft, int direction) {
  int another = direction ^ 1;

  std::queue<std::pair<V, W>> que;
  que.emplace(v_to, d_ft);
  while (!que.empty()) {
    V v;
    W d;
    std::tie(v, d) = que.front();
    que.pop();
    if (query_distance(g, v_from, v, direction) <= d) continue;
    idx[another][v].update(v_from, d);
    for (const auto &w : adj[direction][v]) que.emplace(w, d + 1);
  }
}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(
    const G &g, V v_from, const E &e) {
  V v_to = to(e);
  assert(v_from < g.num_vertices() && v_to < g.num_vertices());
  if (query_distance(g, v_from, v_to) <= 1) return;
  adj[0][v_from].push_back(v_to);
  adj[1][v_to].push_back(v_from);

  std::vector<std::pair<std::pair<V, W>, int>> tmp;
  for (const auto &p : idx[1][v_from].spt) tmp.push_back({p, 0});
  for (const auto &p : idx[0][v_to].spt) tmp.push_back({p, 1});

  std::sort(tmp.begin(), tmp.end());

  for (const auto &q : tmp) {
    auto p = q.first;
    if (q.second == 0) resume_pbfs(g, p.first, v_to, p.second + 1, 0);
    if (q.second == 1) resume_pbfs(g, p.first, v_from, p.second + 1, 1);
  }
}
}  // namespace agl
