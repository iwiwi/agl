#pragma once

#include "graph/graph.h"
#include "heap/heap.h"

namespace agl {
template<typename GraphT = G>
std::vector<typename GraphT::W> single_source_distance(const GraphT &g, V s, D d = kFwd) {
  using W = typename GraphT::W;
  auto h = make_dijkstra_heap(g);
  h.decrease(s, 0);

  while (!h.empty()) {
    V v = h.top_vertex();
    W w = h.top_weight();
    h.pop();
    for (const auto &e : g.edges(v, d)) {
      h.decrease(to(e), w + weight(e));
    }
  }

  return h.weights();
}

template<typename GraphT = G>
std::vector<std::vector<typename GraphT::W>> all_pairs_distance(const GraphT &g, D d = kFwd) {
  using W = typename GraphT::W;
  std::vector<std::vector<W>> dm;

  for (const auto& s : g.vertices()) {
    dm.emplace_back(single_source_distance(g, s, d));
  }

  return dm;
}

std::vector<std::pair<W, double>> single_source_distance_with_num_paths(const G &g, V s, D d = kFwd);
}  // namespace agl
