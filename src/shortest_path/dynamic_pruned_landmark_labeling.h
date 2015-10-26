#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"

namespace agl {
template<size_t kNumBitParallelRoots = 16>
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

  virtual void remove_vertices(const G&, V old_num_vertices) override {
    assert(false);
  }
};

// 確認：テンプレートクラスだから，.cc じゃなくてここに実装を書きます．
template<size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::construct(const G &g) {
  puts("CONSTRUCT ME!");
}

template<size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(const G &g, V v_from, V v_to) {
  puts("hoge");
  return 0;
}

template<size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(const G &g, V v_from, const E &e) {
  puts("new edge");
}
}  // namespace agl
