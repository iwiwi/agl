#pragma once

#include "dijkstra.h"
#include "connectivity/strongly_connected_components.h"

namespace agl {
template<typename GraphT = G>
typename GraphT::W diameter(const GraphT& g, const int kNumDoubleSweep = 10) {
  using W = typename GraphT::W;
  constexpr W kInfW = infinity_weight<W>();
  auto cmp = [](W a, W b){return (a != kInfW ? a : -1) < (b != kInfW ? b : -1);};

  const V num_v = g.num_vertices();
  const auto scc = strongly_connected_components(g);

  // Order vertices
  std::vector<std::tuple<V, double, V>> order;
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
    order.emplace_back(scc[v], -(double)in * out, v);
  }
  std::sort(order.begin(), order.end());

  // Compute the diameter lower bound by the double sweep algorithm
  W diameter = 0;
  for (int i = 0; i < kNumDoubleSweep; ++i) {
    V start = agl::random(num_v);
    auto dist = single_source_distance(g, start);
    start = max_element(dist.begin(), dist.end(), cmp) - dist.begin();
    dist = single_source_distance(g, start, kBwd);
    diameter = std::max(diameter, *max_element(dist.begin(), dist.end(), cmp));
  }

  // Examine every vertex
  std::vector<W> ecc(num_v, kInfW);
  for (V i : make_irange(num_v)) {
    V u = std::get<2>(order[i]);
    if (is_le(ecc[u], diameter)) continue;

    // Refine the eccentricity upper bound
    std::vector<std::pair<V, W>> neighbors;
    for (const auto& e : g.edges(u)) {
      if (ecc[to(e)] != kInfW) {
        neighbors.emplace_back(scc[to(e)], ecc[to(e)] + weight(e));
      } else {
        neighbors.emplace_back(scc[to(e)], kInfW);
      }
    }
    sort(neighbors.begin(), neighbors.end());

    W ub = 0;
    for (size_t j = 0; j < neighbors.size(); ) {
      V component = neighbors[j].first;
      W lb = kInfW;
      for (; j < neighbors.size(); ++j) {
        if (neighbors[j].first != component) break;
        lb = std::min(lb, neighbors[j].second);
      }
      ub = std::max(ub, lb);
      if (is_lt(diameter, ub)) break;
    }

    if (is_le(ub, diameter)) {
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
