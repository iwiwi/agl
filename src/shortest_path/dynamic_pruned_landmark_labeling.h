#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"

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
    std::vector<V> spt_v;
    std::vector<W> spt_d;
  };

  std::vector<index_t> idx[2];

  std::vector<std::vector<V>> adj[2];
  std::vector<V> inv;
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
      double deg = adj[0][i].size() + adj[1][i].size() +
                   (double)agl::random(num_v) / num_v;
      decreasing_order[i] = {deg, i};
    }
    std::sort(decreasing_order.rbegin(), decreasing_order.rend());
    for (int i = 0; i < num_v; ++i) {
      inv[i] = decreasing_order[i].second;
    }
  }

  // Pruned labelling
  for (int i = 0; i < 2; ++i) {
  }
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to) {
  if (v_from >= g.num_vertices() || v_to >= g.num_vertices())
    return v_from == v_to ? 0 : W_INF;

  W d = W_INF;
  const index_t &idx_from = idx[0][v_from];
  const index_t &idx_to = idx[1][v_to];

  for (int i1 = 0, i2 = 0;
       i1 < idx_from.spt_v.size() && i2 < idx_to.spt_v.size();) {
    V v1 = idx_from.spt_v[i1];
    V v2 = idx_from.spt_v[i2];
    if (v1 == v2) {
      W td = idx_from.spt_d[i1] + idx_to.spt_d[i2];
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
