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

 private:
  struct index_t {
    uint32_t *spt_v;
    uint8_t *spt_d;
    uint32_t spt_l;
  } __attribute__((aligned(64)));  // Aligned for cache lines

  index_t *index_;
  std::vector<std::vector<V>> adj_[2];
};

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::construct(
    const G &g) {
  size_t num_edge = g.num_edges();
  size_t num_v = g.num_vertices();
  unweighted_edge_list es = g.edge_list();

  adj_[0].resize(num_v);
  adj_[1].resize(num_v);
  for (size_t i = 0; i < num_edge; ++i) {
    V from = es[i].first, to = es[i].second;
    adj_[0][from].push_back(to);
    adj_[1][to].push_back(from);
  }

  // pruned bfs
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to) {}

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(
    const G &g, V v_from, const E &e) {}
}  // namespace agl
