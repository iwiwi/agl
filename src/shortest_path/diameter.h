#pragma once

#include "dijkstra.h"
#include "connectivity/strongly_connected_components.h"

namespace agl {
template<typename GraphT = G, int kNumDoubleSweep = 10>
typename GraphT::W diameter(const GraphT& g) {
  using W = typename GraphT::W;
  constexpr W INF = std::numeric_limits<W>::max();
  auto cmp = [](W a, W b){return (a != INF ? a : -1) < (b != INF ? b : -1);};

  const V num_v = g.num_vertices();
  const auto scc = strongly_connected_components(g);

  // Order vertices
  std::vector<std::pair<int64_t, V>> order(num_v);
  for (V v : make_irange(num_v)) {
    V in = 0, out = 0;
    for (const auto& e : g.edges(v, kBwd)) {
      if (scc[to(e)] == scc[v]) ++in;
    }
    for (const auto& e : g.edges(v)) {
      if (scc[to(e)] == scc[v]) ++out;
    }
    // SCC : reverse topological order
    // inside an SCC : decreasing order of the product of the indegree and outdegree for vertices in the same SCC
    order[v] = std::make_pair(((int64_t)scc[v] << 32) - in * out, v);
  }
  std::sort(order.begin(), order.end());

  // Compute the diameter lower bound by the double sweep algorithm
  W diameter = 0;
  for (int i = 0; i < kNumDoubleSweep; ++i) {
    V start = agl::random(num_v);
    auto dist = single_source_distance(g, start);
    start = max_element(dist.begin(), dist.end(), cmp) - dist.begin();
    dist = single_source_distance(g, start);
    diameter = std::max(diameter, *max_element(dist.begin(), dist.end(), cmp));
  }

  // Examine every vertex
  std::vector<W> ecc(num_v, INF);
  for (V i : make_irange(num_v)) {
    V u = order[i].second;
    if (ecc[u] <= diameter) continue;

    // Refine the eccentricity upper bound
    W ub = 0;
    std::vector<std::pair<V, W>> neighbors;
    for (const auto& e : g.edges(u)) {
      if (ecc[to(e)] != INF) {
        neighbors.emplace_back(scc[to(e)], ecc[to(e)] + weight(e));
      }
    }
    sort(neighbors.begin(), neighbors.end());

    for (size_t j = 0; j < neighbors.size(); ) {
      V component = neighbors[j].first;
      W lb = INF;
      for (; j < neighbors.size(); ++j) {
        if (neighbors[j].first != component) break;
        lb = std::min(lb, neighbors[j].second);
      }
      ub = std::max(ub, lb);
      if (ub > diameter) break;
    }

    if (ub <= diameter) {
      ecc[u] = ub;
      continue;
    }

    // Update bounds
    auto dist = single_source_distance(g, u);
    ecc[u] = *max_element(dist.begin(), dist.end(), cmp);
    diameter = std::max(diameter, ecc[u]);

    auto h = make_dijkstra_heap(g);
    h.decrease(u, 0);
    while (!h.empty()) {
      V v = h.top_vertex();
      W w = h.top_weight();
      h.pop();
      ecc[v] = std::min(ecc[v], w + ecc[u]);
      for (const auto &e : g.edges(v, kBwd)) {
        if (scc[to(e)] == scc[u]) {
          h.decrease(to(e), w + weight(e));
        }
      }
    }
  }

  return diameter;
}
}  // namespace agl
