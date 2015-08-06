#pragma once
#include "radix_heap.h"
#include "graph/graph.h"

namespace agl {
template<typename GraphT = G>
class dijkstra_heap {
public:
  using V = typename GraphT::V;
  using W = typename GraphT::W;

  explicit dijkstra_heap(V num_vertices = 0) : ws_(num_vertices, kInfW) {}

  inline bool decrease(V v, W w) {
    if (is_le(ws_[v], w)) return false;
    ws_[v] = w;
    h_.push(w, v);
    canonicalize();
    return true;
  }

  inline void pop() {
    vs_.emplace_back(top_vertex());
    h_.pop();
    canonicalize();
  }

  void clear() {
    for (auto v : vs_) ws_[v] = kInfW;
    vs_.clear();
    while (!h_.empty()) {
      ws_[top_vertex()] = kInfW;
      h_.pop();
    }
    h_.clear();
  }

  inline bool empty() const {
    return h_.empty();
  }

  inline V top_vertex()  {
    return h_.top_value();
  }

  inline W top_weight() {
    return h_.top_key();
  }

  inline const std::vector<W> &weights() const {
    return ws_;
  }

private:
  using internal_heap = radix_heap::pair_radix_heap<W, V>;

  std::vector<W> ws_;
  internal_heap h_;
  std::vector<V> vs_;

  void canonicalize() {
    while (!h_.empty()) {
      if (ws_[h_.top_value()] == h_.top_key()) break;
      h_.pop();
    }
  }
};

template<typename GraphT>
dijkstra_heap<GraphT> make_dijkstra_heap(const GraphT &g) {
  return dijkstra_heap<GraphT>(g.num_vertices());
}
}
