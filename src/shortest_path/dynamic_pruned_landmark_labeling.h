#pragma once
#include "graph/graph.h"
#include "graph/graph_index_interface.h"
#include <malloc.h>

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
  const uint8_t INF8 = 100;  // For unreachable pairs
  const uint32_t INF32 = std::numeric_limits<int32_t>::max();  // For sentinel

  struct index_t {
    uint32_t *spt_vertex;
    uint8_t *spt_dist;
    uint32_t spt_l;
  } __attribute__((aligned(64)));  // Aligned for cache lines

  index_t *index_, *index_inv_;
  std::vector<std::vector<V>> adj_, adj_inv_;
  V num_vertices_;
  E num_edges_;
};

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::construct(
    const G &g) {
  E &num_e = num_edges_;
  V &num_v = num_vertices_;
  num_e = g.num_edges();
  num_v = g.num_vertices();
  unweighted_edge_list es = g.edge_list();

  adj_.resize(num_v);
  adj_inv_.resize(num_v);
  for (size_t i = 0; i < num_v; ++i) {
    V from = es[i].first, to = es[i].second;
    adj_[from].push_back(to);
    adj_inv_[to].push_back(from);
  }

  // Initialize Index

  index_ = (index_t *)memalign(64, num_v * sizeof(index_t));
  index_inv_ = (index_t *)memalign(64, num_v * sizeof(index_t));
  if (index_ == NULL || index_inv_ == NULL) {
    num_v = 0;
    return;
  }
  for (V v = 0; v < num_v; ++v) {
    index_[v].spt_vertex = NULL;
    index_[v].spt_dist = NULL;
    index_[v].spt_l = 0;
    index_inv_[v].spt_vertex = NULL;
    index_inv_[v].spt_dist = NULL;
    index_inv_[v].spt_l = 0;
  }

  // Pruned bfs
  {
    std::vector<W> tmp_dist(num_v);
    std::vector<V> queue(num_v);
    for (V v = 0; v < num_v; ++v) {
      int queue_tail = 0, queue_head = 0;
      queue[queue_tail++] = v;

      fill(tmp_dist.begin(), tmp_dist.end(), INF8);
      tmp_dist[v] = 0;
      // L<-L
      while (queue_head < queue_tail) {
        V u = queue[queue_head++];
        
      }
    }
  }
}

template <size_t kNumBitParallelRoots>
W dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::query_distance(
    const G &g, V v_from, V v_to) {
  if (v_from == v_to) return 0;
  if (v_from >= num_vertices_ || v_to >= num_vertices_) return kInfW;

  const index_t &idx_from = index_[v_from];
  const index_t &idx_to = index_inv_[v_to];
  W dist = INF8;

  for (int i1 = 0, i2 = 0;;) {
    V v1 = idx_from.spt_vertex[i1], v2 = idx_to.spt_vertex[i2];
    if (v1 == v2) {
      if (v1 == INF32) break;  // Sentinel
      W tmp_dist = idx_from.spt_dist[i1] + idx_to.spt_dist[i2];
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

template <size_t kNumBitParallelRoots>
void dynamic_pruned_landmark_labeling<kNumBitParallelRoots>::add_edge(
    const G &g, V v_from, const E &e) {}
}  // namespace agl
